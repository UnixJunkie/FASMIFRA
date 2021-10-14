(* Copyright (C) 2021, Francois Berenger
   Tsuda Laboratory,
   Tokyo University,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Ultra fast (Deep)SMILES molecular fragments assembler *)

open Printf

module A = Array
module Buff = Buffer
module CLI = Minicli.CLI
module Fn = Filename
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log

type input_smi_token = Open_paren
                     | Close_paren
                     | Cut_bond of int * int
                     | Bracket_atom of string
                     | Ring_closure of int
                     | Rest of string

let parens_regexp = Str.regexp "[()]"
let cut_bond_regexp = Str.regexp "\\[[0-9]+\\*\\]\\[[0-9]+\\*\\]"
(* cut bond with a parenthesis inside must be rewritten *)
let paren_cut_bond_regexp = Str.regexp "\\[[0-9]+\\*\\](\\[[0-9]+\\*\\]"
let bracket_atom_regexp = Str.regexp "\\[[^]]+\\]"
let double_digits_regexp = Str.regexp "%[0-9][0-9]"
let single_digit_regexp = Str.regexp "[0-9]"

let renumber_ring_closure ht (smi_depth: int) (ring_closure: int) =
  try Ht.find ht (smi_depth, ring_closure)
  with Not_found ->
    let res = Ht.length ht in
    Ht.add ht (smi_depth, ring_closure) res;
    res

let parse_cut_bond s =
  try Scanf.sscanf s "[%d*][%d*]" (fun i j -> Cut_bond (i, j))
  with exn ->
    let () = Log.fatal "Fasmifra.parse_cut_bond: cannot parse '%s'" s in
    raise exn

let parse_paren_cut_bond s =
  try Scanf.sscanf s "[%d*]([%d*]" (fun i j -> (i, j))
  with exn ->
    let () = Log.fatal "Fasmifra.parse_paren_cut_bond: cannot parse '%s'" s in
    raise exn

let string_of_ring_closure c =
  if c > 99 then (* forbid triple digits ring closure *)
    (Log.fatal "Fasmifra.string_of_ring_closure: %d > 99" c;
     exit 1)
  else if c > 9 then
    (* a double digit ring closure in SMILES must be prefixed by % *)
    sprintf "%%%d" c
  else
    sprintf "%d" c

let fprintf_ring_closure out c =
  if c > 99 then (* forbid triple digits ring closure *)
    (Log.fatal "Fasmifra.fprintf_ring_closure: %d > 99" c;
     exit 1)
  else if c > 9 then
    fprintf out "%%%d" c
  else
    fprintf out "%d" c

let bprintf_ring_closure buff c =
  if c > 99 then (* forbid triple digits ring closure *)
    (Log.fatal "Fasmifra.bprintf_ring_closure: %d > 99" c;
     exit 1)
  else if c > 9 then
    bprintf buff "%%%d" c
  else
    bprintf buff "%d" c

let string_of_smi_token = function
  | Open_paren -> "("
  | Close_paren -> ")"
  | Cut_bond (i, j) ->
    (* only OK in debugging output: those should disappear
     * during fragments assembly *)
    sprintf "[%d*][%d*]" i j
  | Ring_closure rc -> string_of_ring_closure rc
  | Bracket_atom x
  | Rest x -> x

(* for performance *)
let fprintf_smi_token out = function
  | Open_paren -> output_char out '('
  | Close_paren -> output_char out ')'
  | Cut_bond _ -> assert(false) (* should have been deleted earlier *)
  | Ring_closure rc -> fprintf_ring_closure out rc
  | Bracket_atom x
  | Rest x -> output_string out x

(* for the parallel version, if we let the workers do the string conversion *)
(* it might be faster to let the muxer to do it though, TO TEST *)
let bprintf_smi_token buff = function
  | Open_paren -> Buff.add_char buff '('
  | Close_paren -> Buff.add_char buff ')'
  | Cut_bond _ -> assert(false) (* should have been deleted earlier *)
  | Ring_closure rc -> bprintf_ring_closure buff rc
  | Bracket_atom x
  | Rest x -> Buff.add_string buff x

let string_of_tokens tokens =
  String.concat "" (L.map string_of_smi_token tokens)

let fprintf_tokens out tokens =
  L.iter (fprintf_smi_token out) tokens

let sprintf_tokens tokens =
  let buff = Buff.create 256 in
  L.iter (bprintf_smi_token buff) tokens;
  Buff.contents buff

let parse_double_digit_ring_closure s =
  try Scanf.sscanf s "%%%d" (fun x -> assert(x >= 10); x)
  with exn ->
    (Log.fatal "Fasmifra.parse_double_digit_ring_closure: cannot parse: '%s'" s;
     raise exn)

let parse_single_digit_ring_closure s =
  try Scanf.sscanf s "%d" (fun x -> assert(x <= 9); x)
  with exn ->
    (Log.fatal "Fasmifra.parse_single_digit_ring_closure: cannot parse: '%s'" s;
     raise exn)

let tokenize_full (s: string): input_smi_token list =
  let res =
    L.map Str.(function
        | Delim "(" -> [[[[Open_paren]]]]
        | Delim ")" -> [[[[Close_paren]]]]
        | Delim _ -> assert(false) (* parens_regexp would be wrong then *)
        | Text a ->
          L.map Str.(function
              | Delim cut_bond_str -> [[[parse_cut_bond cut_bond_str]]]
              | Text b ->
                L.map Str.(function
                    | Delim bracket_atom_str -> [[Bracket_atom bracket_atom_str]]
                    | Text c ->
                      L.map Str.(function
                          | Delim double_digit_rc ->
                            [Ring_closure (parse_double_digit_ring_closure double_digit_rc)]
                          | Text d ->
                            L.map Str.(function
                                | Delim single_digit_rc ->
                                  Ring_closure (parse_single_digit_ring_closure single_digit_rc)
                                | Text e -> Rest e
                              )
                              (Str.bounded_full_split single_digit_regexp d 1024)
                        )
                        (Str.bounded_full_split double_digits_regexp c 1024)
                  )
                  (Str.bounded_full_split bracket_atom_regexp b 1024)
            )
            (Str.bounded_full_split cut_bond_regexp a 1024)
      )
      (* WARNING: bug if more than 1024 tokens in the string ! *)
      (Str.bounded_full_split parens_regexp s 1024) in
  L.flatten (L.flatten (L.flatten (L.flatten res)))

(* like BatList.takedrop but considering parenthesis nesting depth
   take as long as the current depth is >= to the initial one *)
let take_drop_depth depth tokens =
  (* Log.warn "take_drop_depth: %d" depth; *)
  let rec loop d acc = function
    | [] -> (L.rev acc, [])
    | tok :: toks ->
      match tok with
      | Open_paren -> loop (d + 1) (tok :: acc) toks
      | Close_paren ->
        if d = depth then
          (L.rev acc, tok :: toks)
        else
          loop (d - 1) (tok :: acc) toks
      | _ -> loop d (tok :: acc) toks
  in
  if depth = 0 then (tokens, [])
  else loop depth [] tokens

(* return the current fragment, plus a list of token lists to fragment
 * further *)
let fragment_once (tokens: input_smi_token list):
  (input_smi_token list) * (input_smi_token list list) =
  (* Log.warn "fragment_once: %s" (string_of_tokens tokens); *)
  let rec loop depth (seed, branches) = function
    | [] -> (L.rev seed, branches)
    | tok :: toks ->
      match tok with
      | Open_paren -> loop (depth + 1) (tok :: seed, branches) toks
      | Close_paren -> loop (depth - 1) (tok :: seed, branches) toks
      | Cut_bond (i, j) ->
        let same_depth, rest = take_drop_depth depth toks in
        (* Log.info "same: %s" (string_of_tokens same_depth);
         * Log.info "rest: %s" (string_of_tokens rest); *)
        let branch = Cut_bond (j, i) :: same_depth in
        loop depth (tok :: seed, branch :: branches) rest
      | _ -> loop depth (tok :: seed, branches) toks
  in
  match tokens with
  | Cut_bond (i, j) :: rest ->
    (* in a branch, the leading cut bond must be discarded to avoid
       an infinite loop;
       flipping (i, j) to (j, i) is required *)
    loop 0 ([Cut_bond (j, i)], []) rest
  | _ -> loop 0 ([], []) tokens

let rec fragment seeds ht tokens =
  (* Log.warn "fragment_____: %s" (string_of_tokens tokens); *)
  let seed, branches = fragment_once tokens in
  (* Log.info "seed: %s" (string_of_tokens seed);
   * L.iter (fun branch ->
   *     Log.info "branch: %s" (string_of_tokens branch)
   *   ) branches; *)
  begin
    match seed with
    | [] -> assert(false)
    | x :: xs ->
      match x with
      | Cut_bond (i, j) ->
        (* a branch starts with a cut bond *)
        let prev_branches = Ht.find_default ht (i, j) [] in
        (* typed cut bond discarded here *)
        Ht.replace ht (i, j) (xs :: prev_branches)
      | _ ->
        (* a seed does not start with a cut bond *)
        seeds := seed :: !seeds
  end;
  L.iter (fragment seeds ht) branches

(* handle rare case: put the parenthesis before the first atom.
   e.g: "[1*]([2*]" --> "([1*][2*]". Not sure this is still needed. *)
let rewrite_paren_cut_bond_smiles s =
  String.concat ""
    (L.map Str.(function
         | Delim paren_cut_bond ->
           let (i, j) = parse_paren_cut_bond paren_cut_bond in
           sprintf "([%d*][%d*]" i j
         | Text t -> t
       ) (Str.bounded_full_split paren_cut_bond_regexp s 1024)
    )

let index_fragments named_smiles =
  let n = L.length named_smiles in
  let seed_frags = ref [] in
  let frags_ht = Ht.create n in
  L.iter (fun (smiles, _name) ->
      let rewritten = rewrite_paren_cut_bond_smiles smiles in
      let tokens = tokenize_full rewritten in
      fragment seed_frags frags_ht tokens
    ) named_smiles;
  let seeds = A.of_list !seed_frags in
  let ht = Ht.map (fun _ij branches -> A.of_list branches) frags_ht in
  (seeds, ht)

let rev_renumber_ring_closures ht i tokens =
  let res =
    L.rev_map (function
        | Ring_closure j -> Ring_closure (renumber_ring_closure ht !i j)
        | tok -> tok
      ) tokens in
  incr i;
  res

let array_rand_elt rng a =
  let n = A.length a in
  let i = BatRandom.State.int rng n in
  Array.unsafe_get a i

let assemble_smiles_fragments rng seeds branches =
  let frag_count = ref 0 in
  let ht = Ht.create 97 in
  let seed_frag =
    let chosen = array_rand_elt rng seeds in
    (* Log.error "chosen seed: %s" (string_of_tokens chosen); *)
    L.rev (rev_renumber_ring_closures ht frag_count chosen) in
  (* Log.error "renumbered seed: %s" (string_of_tokens seed_frag); *)
  let rec loop acc tokens = match tokens with
    | [] -> L.rev acc
    | x :: xs ->
      match x with
      | Cut_bond (i, j) ->
        let possible_branches = Ht.find branches (i, j) in
        let branch =
          let chosen = array_rand_elt rng possible_branches in
          (* Log.error "chosen branch: %s" (string_of_tokens chosen); *)
          rev_renumber_ring_closures ht frag_count chosen in
        (* Log.error "renumbered branch: %s" (string_of_tokens (L.rev branch)); *)
        (* typed cut bond discarded here *)
        loop acc (L.rev_append branch xs)
      | _  -> loop (x :: acc) xs
  in
  loop [] seed_frag

(* almost copy/paste of assemble_smiles_fragments *)
let assemble_deepsmiles_fragments rng seeds branches =
  let rec loop acc tokens = match tokens with
    | [] -> L.rev acc
    | x :: xs ->
      match x with
      | Cut_bond (i, j) ->
        let possible_branches = Ht.find branches (i, j) in
        let branch = array_rand_elt rng possible_branches in
        (* typed cut bond discarded here *)
        loop acc (L.append branch xs)
      | _  -> loop (x :: acc) xs
  in
  let seed_frag = array_rand_elt rng seeds in
  loop [] seed_frag

type rng_style = Performance of BatRandom.State.t
               | Repeatable of BatRandom.State.t

(* [Random.int bound]: bound must be greater than 0 and less than 2{^30}. *)
let rand_int_bound = (BatInt.pow 2 30) - 1

let get_rng = function
  | Performance _ -> failwith "Fasmifra.get_rng: don't call this if you want performance"
  | Repeatable rng ->
    let seed = BatRandom.State.int rng rand_int_bound in
    BatRandom.State.make [|seed|]

(* tell how much the worker must generate *)
let demux total current csize () =
  let to_generate = total - !current in
  if to_generate = 0 then
    raise Parany.End_of_input
  else
    let curr_csize = min csize to_generate in
    current := !current + curr_csize;
    curr_csize

(* just to get the type right *)
let worker_rng = ref (BatRandom.State.make_self_init ())

(* setup this worker's RNG to be (very probably) uniq *)
let init rng rank =
  let seed = ref 0 in
  for _i = 0 to rank do
    seed := BatRandom.State.int rng rand_int_bound
  done;
  Log.info "rank: %d seed: %d" rank !seed; (* visual check they differ *)
  worker_rng := BatRandom.State.make [|!seed|]

let work f seed_fragments frags_ht n =
  Array.init n (fun _i ->
      sprintf_tokens (f !worker_rng seed_fragments frags_ht)
    )

let mux out count tokens_str_a =
  A.iter (fun token_str ->
      fprintf out "%s\tgenmol_%d\n" token_str !count;
      incr count
    ) tokens_str_a

let parse_SMILES_line line =
  try Scanf.sscanf line "%s@\t%s" (fun smiles name -> (smiles, name))
  with exn -> (Log.fatal "parse_SMILES_line: malformed line: '%s'" line;
               raise exn)

(* marshal x to file *)
let save fn x =
  LO.with_out_file fn (fun out ->
      Marshal.to_channel out x [Marshal.No_sharing]
    )

(* unmarshal x from file *)
let restore fn =
  LO.with_in_file fn Marshal.from_channel

let cache_indexed_fragments force frags_fn seed_frags_frags_ht_pair =
  let cache_fn = frags_fn ^ ".bin_cache" in
  if not (Sys.file_exists cache_fn) || force then
    let () =
      if force then
        Log.warn "overwriting indexed fragments cache: %s" cache_fn
      else
        Log.info "creating indexed fragments cache: %s" cache_fn in
    save cache_fn seed_frags_frags_ht_pair
  else
    Log.warn "cache file already exists (use -f to overwrite): %s" cache_fn

let load_indexed_fragments force frags_fn =
  let cache_fn = frags_fn ^ ".bin_cache" in
  if (not force) && Sys.file_exists cache_fn then
    let () = Log.info "reading indexed fragments from cache: %s" cache_fn in
    restore cache_fn
  else
    let input_frags = LO.map frags_fn parse_SMILES_line in
    index_fragments input_frags

let main () =
  let start = Unix.gettimeofday () in
  (* Logger ---------------------------------------------------------------- *)
  Log.color_on ();
  Log.(set_log_level INFO);
  Log.(set_prefix_builder short_prefix_builder);
  (* CLI ------------------------------------------------------------------- *)
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              -i <filename>: smiles fragments input file\n  \
              -o <filenams>: output file\n  \
              -n <int>: how many molecules to generate\n  \
              [-f]: overwrite existing indexed fragments cache file\n  \
              [-np <int>]: max number of processes (default=1)\n  \
              [-c <int>]: chunk size (for -np; default=1)\n  \
              [--seed <int>]: RNG seed\n  \
              [--deep-smiles]: input/output molecules in DeepSMILES\n  \
              no-rings format\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let n = CLI.get_int ["-n"] args in
  (* to make things truly repeatable, this rng is in fact a stream
   * of RNG seeds *)
  let rng_style = match CLI.get_int_opt ["--seed"] args with
    | None -> Performance (BatRandom.State.make_self_init ())
    | Some seed -> Repeatable (BatRandom.State.make [|seed|]) in
  let input_frags_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let csize = CLI.get_int_def ["-c"] args 1 in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let force = CLI.get_set_bool ["-f"] args in
  let assemble =
    if CLI.get_set_bool ["--deep-smiles"] args then
      assemble_deepsmiles_fragments
    else
      assemble_smiles_fragments in
  CLI.finalize (); (* ------------ CLI parsing ---------------- *)
  Log.info "indexing fragments";
  let seed_fragments, frags_ht =
    let res = load_indexed_fragments force input_frags_fn in
    cache_indexed_fragments force input_frags_fn res;
    res in
  Log.info "seed_frags: %d; attach_types: %d"
    (A.length seed_fragments) (Ht.length frags_ht);
  LO.with_out_file output_fn (fun out ->
      if nprocs <= 1 then
        match rng_style with
        | Performance rng ->
          for i = 1 to n do
            let tokens = assemble rng seed_fragments frags_ht in
            fprintf_tokens out tokens;
            fprintf out "\tgenmol_%d\n" i
          done
        | Repeatable _ ->
          for i = 1 to n do
            let rng = get_rng rng_style in
            let tokens = assemble rng seed_fragments frags_ht in
            fprintf_tokens out tokens;
            fprintf out "\tgenmol_%d\n" i
          done
      else (* nprocs > 1 then *)
        match rng_style with
        | Repeatable _ ->
          failwith "only Performance rng supported w/ nprocs > 1"
        | Performance rng ->
          Parany.run
            ~init:(init rng)
            ~csize:1 (* fixed; keep 1 here *)
            nprocs
            ~demux:(demux n (ref 0) csize)
            ~work:(work assemble seed_fragments frags_ht)
            ~mux:(mux out (ref 0))
    );
  let stop = Unix.gettimeofday () in
  let dt = stop -. start in
  Log.info "rate: %.2f molecule/s" ((float n) /. dt)

let () = main ()
