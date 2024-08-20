(* Copyright (C) 2024, Francois Berenger
   Tsuda Laboratory,
   The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Ultra fast (Deep)SMILES molecular fragments assembler *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fn = Filename
module Ht = BatHashtbl
module ISet = BatSet.Int
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

type input_smi_token = Open_paren
                     | Close_paren
                     | Cut_bond of int * int
                     | Bracket_atom of string
                     | Ring_closure of int
                     | Rest of string

let parens_regexp = Str.regexp "[()]"
let cut_bond_regexp = Str.regexp "\\[\\*:[0-9]+\\]\\[\\*:[0-9]+\\]"
(* cut bond with a parenthesis inside must be rewritten *)
let paren_cut_bond_regexp = Str.regexp "\\[\\*:[0-9]+\\](\\[\\*:[0-9]+\\]"
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
  try Scanf.sscanf s "[*:%d][*:%d]" (fun i j -> Cut_bond (i, j))
  with exn ->
    let () = Log.fatal "Fasmifra.parse_cut_bond: cannot parse '%s'" s in
    raise exn

let parse_paren_cut_bond s =
  try Scanf.sscanf s "[*:%d]([*:%d]" (fun i j -> (i, j))
  with exn ->
    let () = Log.fatal "Fasmifra.parse_paren_cut_bond: cannot parse '%s'" s in
    raise exn

exception Too_many_rings

let string_of_ring_closure c =
  if c > 99 then (* forbid triple digits ring closure *)
    (Log.fatal "Fasmifra.string_of_ring_closure: %d > 99" c;
     raise Too_many_rings)
  else if c > 9 then
    (* a double digit ring closure in SMILES must be prefixed by % *)
    sprintf "%%%d" c
  else
    sprintf "%d" c

let fprintf_ring_closure out c =
  if c > 99 then (* forbid triple digits ring closure *)
    (Log.fatal "Fasmifra.fprintf_ring_closure: %d > 99" c;
     raise Too_many_rings)
  else if c > 9 then
    fprintf out "%%%d" c
  else
    fprintf out "%d" c

(* for debugging *)
let string_of_smi_token = function
  | Open_paren -> "("
  | Close_paren -> ")"
  | Cut_bond (i, j) ->
    (* only OK in debugging output: those should disappear
       during fragments assembly *)
    sprintf "[*:%d][*:%d]" i j
  | Ring_closure rc -> string_of_ring_closure rc
  | Bracket_atom x
  | Rest x -> x

(* for performance *)
let fprintf_smi_token out = function
  | Open_paren -> output_char out '('
  | Close_paren -> output_char out ')'
  | Cut_bond (i, j) ->
    (* in some use cases, cut bonds are supposed to have been deleted by
       assemble_smiles_fragments *)
    fprintf out "[*:%d][*:%d]" i j
  | Ring_closure rc -> fprintf_ring_closure out rc
  | Bracket_atom x
  | Rest x -> output_string out x

let string_of_tokens tokens =
  String.concat "" (L.map string_of_smi_token tokens)

let fprintf_tokens out tokens =
  L.iter (fprintf_smi_token out) tokens

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
              | Delim cut_bond_str ->
                [[[parse_cut_bond cut_bond_str]]]
              | Text b ->
                L.map Str.(function
                    | Delim bracket_atom_str ->
                      [[Bracket_atom bracket_atom_str]]
                    | Text c ->
                      L.map Str.(function
                          | Delim double_digit_rc ->
                            [Ring_closure (parse_double_digit_ring_closure
                                             double_digit_rc)]
                          | Text d ->
                            L.map Str.(function
                                | Delim single_digit_rc ->
                                  Ring_closure (parse_single_digit_ring_closure
                                                  single_digit_rc)
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
   e.g: "[*:1]([*:2]" --> "([*:1][*:2]". Not sure this is still needed. *)
let rewrite_paren_cut_bond_smiles s =
  String.concat ""
    (L.map Str.(function
         | Delim paren_cut_bond ->
           let (i, j) = parse_paren_cut_bond paren_cut_bond in
           sprintf "([*:%d][*:%d]" i j
         | Text t -> t
       ) (Str.bounded_full_split paren_cut_bond_regexp s 1024)
    )

(* dump all fragments to opened file.
   REMARK: each fragment is a valid SMILES if cut bonds are not erased. *)
let dump_seed_fragments
    (out: out_channel) (frags: input_smi_token list array): unit =
  A.iter (fun tokens ->
      fprintf_tokens out tokens;
      output_char out '\n' (* terminate fragment *)
    ) frags

(* almost like dump_seed_fragments *)
let dump_branch_fragments
    (out: out_channel)
    (cut_bond: input_smi_token)
    (frags: input_smi_token list array): unit =
  A.iter (fun tokens ->
      (* properly outputing a branch fragment requires
         that it is prefixed by the correct cut bond *)
      fprintf_tokens out (cut_bond :: tokens);
      output_char out '\n' (* terminate fragment *)
    ) frags

let index_fragments maybe_out_fn named_smiles =
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
  (match maybe_out_fn with
   | None -> ()
   | Some output_fn ->
     LO.with_out_file output_fn (fun out ->
         dump_seed_fragments out seeds;
         Ht.iter (fun (i, j) branches ->
             dump_branch_fragments out (Cut_bond (i, j)) branches
           ) ht
       );
     Log.info "SMILES fragments written to %s" output_fn;
     exit 0 (* do not try to generate molecules after this *)
  );
  (seeds, ht)

(* unique identifier for each fragment *)
type frag_id = { i_j: int * int; (* cut bond type;
                                    (-1,-1) for seeds by convention *)
                 k: int } (* index in array for that cut bond type *)

let string_of_frag_id x =
  let i, j = x.i_j in
  sprintf "%d-%d:%d" i j x.k

let string_of_frag_ids l =
  (* a molecule is supposed to be made only of a handful
     of fragments, so this is OK *)
  S.concat "," (L.rev_map string_of_frag_id l)

let frag_id_of_string x =
  Scanf.sscanf x "%d-%d:%d" (fun i j k -> { i_j = (i, j); k })

let frag_ids_of_string s =
  L.map frag_id_of_string (S.split_on_char ',' s)

(* a Gaussian distribution *)
type dist = { mu: float;
              sigma: float }

let string_of_dist d =
  (* think about m+/-s, but shorter *)
  sprintf "%g/%g" d.mu d.sigma

let dist_of_string s =
  Scanf.sscanf s "%f/%f" (fun mu sigma -> { mu; sigma })

(* attach a distribution to each seed and each branch fragment
   we should do this only at the first iteration:
   -mu and -sigma are provided; no -ig
   in subsequent iterations: -sigma and -ig are provided *)
let initialize_gaussians seeds_a branches_ht mu sigma =
  let init_dist = { mu; sigma } in
  let res = Ht.create (1 + Ht.length branches_ht) in
  (* gaussians for seed fragments *)
  Ht.add res (-1, -1) (let num_seeds = A.length seeds_a in
                       A.make num_seeds init_dist);
  (* gaussians for branch fragments *)
  Ht.iter (fun i_j branches ->
      let num_frags = A.length branches in
      Ht.add res i_j (A.make num_frags init_dist)
    ) branches_ht;
  res

let rev_renumber_ring_closures ht i tokens =
  let res =
    L.rev_map (function
        | Ring_closure j -> Ring_closure (renumber_ring_closure ht !i j)
        | tok -> tok
      ) tokens in
  incr i;
  res

let almost_one = Float.pred 1.0
(* enforce that the highest float we can draw is <1.0 *)
let _ = assert(almost_one < 1.0)
let pi = 4.0 *. (atan 1.0)
let two_pi = 2.0 *. pi

(* [gauss mu sigma] get one float from the normal distribution
   with mean=mu and stddev=sigma
   a = cos(2*pi*x) * sqrt(-2*log(1-y))
   b = sin(2*pi*x) * sqrt(-2*log(1-y)) (b is ignored below)
   cf. Python's documentation of random.gauss function *)
let gauss rng dist =
  dist.mu +. (dist.sigma *.
              (cos (two_pi *. (Random.State.float rng almost_one)) *.
               sqrt (-2.0 *. log (1.0 -. (Random.State.float rng almost_one)))))

(* formulas (4) and (5) in p. 1160 of
 * "Thompson Sampling-An Efficient Method for Searching Ultralarge Synthesis
 * on Demand Databases"; https://doi.org/10.1021/acs.jcim.3c01790
 * [d_t]: distribution for fragment i at time t
 * [s2]: global variance for all fragments (initial estimate)
 * [x_t]: observation for fragment i at t
 * returns the updated distribution for fragment i at t+1 *)
let update_gaussian d_t s2 x_t =
  let s2_t = d_t.sigma *. d_t.sigma in
  let denom = s2_t +. s2 in
  { mu = ((s2_t *. x_t) +. (s2 *. d_t.mu)) /. denom;
    sigma = (s2_t *. s2) /. denom }

(* for one molecule whose composition is known,
 * (seed and fragments ids), update impacted beliefs *)
let update_gaussians dists_ht s2 score frag_ids =
  L.iter (fun frag_id ->
      let arr = Ht.find dists_ht frag_id.i_j in
      arr.(frag_id.k) <- update_gaussian arr.(frag_id.k) s2 score
    ) frag_ids

let update_many_gaussians s2 ij2dists smi_name_scores =
  L.iter (fun (_smi, mol_name, score) ->
      let frag_ids = frag_ids_of_string mol_name in
      update_gaussians ij2dists s2 score frag_ids
    ) smi_name_scores

(* default fragment sampling policy for training-set distribution matching *)
let uniform_random rng _ij n =
  BatRandom.State.int rng n

(* for maximization *)
let thompson_max ij2dists rng ij _n =
  let dists = Ht.find ij2dists ij in
  let samples = A.map (gauss rng) dists in
  let maxi = A.max samples in
  (* uniform random choice among fragments that scored best *)
  let candidates = ref [] in
  A.iteri (fun i sample ->
      if sample = maxi then
        candidates := i :: !candidates
    ) samples;
  let cands = A.of_list !candidates in
  let i = uniform_random rng (-1,-1) (A.length cands) in
  A.unsafe_get cands i

let assemble_smiles_fragments choose_frag_idx rng seeds branches =
  let frag_count = ref 0 in
  let ht = Ht.create 97 in
  let seed_frag =
    let k = choose_frag_idx rng (-1,-1) (A.length seeds) in
    let chosen = A.unsafe_get seeds k in
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
          let k = choose_frag_idx rng (i, j) (A.length possible_branches) in
          let chosen = A.unsafe_get possible_branches k in
          (* Log.error "chosen branch: %s" (string_of_tokens chosen); *)
          rev_renumber_ring_closures ht frag_count chosen in
        (* Log.error "renumbered branch: %s" (string_of_tokens (L.rev branch)); *)
        (* typed cut bond discarded here *)
        loop acc (L.rev_append branch xs)
      | _  -> loop (x :: acc) xs
  in
  (loop [] seed_frag, [])

(* like assemble_smiles_fragments, but "Preserve Cut Bonds" (PCB).
   To output generated molecules w/ cut bonds preserved; so that
   generated molecules do not need to be fragmented later on.
   Also, keep track of the fragments making the molecule. *)
let assemble_smiles_fragments_PCB choose_frag_idx rng seeds branches =
  (* keep track of the fragments making the molecule *)
  let frag_ids = ref [] in
  let frag_count = ref 0 in
  let ht = Ht.create 97 in
  let seed_frag =
    let k = choose_frag_idx rng (-1,-1) (A.length seeds) in
    let chosen = A.unsafe_get seeds k in
    frag_ids := { i_j = (-1,-1); k } :: !frag_ids;
    L.rev (rev_renumber_ring_closures ht frag_count chosen) in
  let rec loop acc tokens = match tokens with
    | [] -> L.rev acc
    | x :: xs ->
      match x with
      | Cut_bond (i, j) ->
        let possible_branches = Ht.find branches (i, j) in
        let branch =
          let k = choose_frag_idx rng (i, j) (A.length possible_branches) in
          let chosen = A.unsafe_get possible_branches k in
          frag_ids := { i_j = (i,j); k } :: !frag_ids;
          rev_renumber_ring_closures ht frag_count chosen in
        (* preserve cut bond [x] here *)
        loop (x :: acc) (L.rev_append branch xs)
      | _  -> loop (x :: acc) xs
  in
  (loop [] seed_frag, !frag_ids)

(* almost copy/paste of assemble_smiles_fragments *)
let assemble_deepsmiles_fragments choose_frag_idx rng seeds branches =
  let rec loop acc tokens = match tokens with
    | [] -> L.rev acc
    | x :: xs ->
      match x with
      | Cut_bond (i, j) ->
        let possible_branches = Ht.find branches (i, j) in
        let k = choose_frag_idx rng (i, j) (A.length possible_branches) in
        let branch = A.unsafe_get possible_branches k in
        (* typed cut bond discarded here *)
        loop acc (L.append branch xs)
      | _  -> loop (x :: acc) xs
  in
  let k = choose_frag_idx rng (-1,-1) (A.length seeds) in
  let seed_frag = A.unsafe_get seeds k in
  (loop [] seed_frag, [])

type rng_style = Performance
               | Repeatable

let parse_SMILES_line line =
  try Scanf.sscanf line "%s@\t%s" (fun smiles name -> (smiles, name))
  with exn -> (Log.fatal "parse_SMILES_line: malformed line: '%s'" line;
               raise exn)

let parse_score_line line =
  try Scanf.sscanf line "%s@\t%f" (fun name score -> (name, score))
  with exn -> (Log.fatal "parse_score_line: malformed line: '%s'" line;
               raise exn)

(* zip named SMILES with named scores *)
let load_scores (smi_fn: string) (scores_fn: string):
  (string * string * float) list =
  let smi_names = LO.map smi_fn parse_SMILES_line in
  let name_scores = LO.map scores_fn parse_score_line in
  L.map2 (fun (smi, name) (name', score) ->
      assert(name = name');
      (smi, name, score)
    ) smi_names name_scores

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

let load_indexed_fragments maybe_out_fn force frags_fn =
  let cache_fn = frags_fn ^ ".bin_cache" in
  if (not force) && Sys.file_exists cache_fn then
    let () = Log.info "reading indexed fragments from cache: %s" cache_fn in
    restore cache_fn
  else
    let input_frags = LO.map frags_fn parse_SMILES_line in
    index_fragments maybe_out_fn input_frags

(* fragments dictionary header line *)
let dict_header = "#i_j:mean_stddevs"

(* save distributions to file *)
let save_gaussians ht fn =
  LO.with_out_file fn (fun out ->
      fprintf out "%s\n" dict_header;
      Ht.iter (fun i_j arr ->
          let i, j = i_j in
          fprintf out "%d-%d:" i j;
          A.iteri (fun k dist ->
              if k > 0 then
                fprintf out ",%g/%g" dist.mu dist.sigma
              else
                fprintf out "%g/%g" dist.mu dist.sigma
            ) arr
        ) ht
    )

let load_gaussians fn =
  let lines = LO.lines_of_file fn in
  match lines with
  | [] -> assert false
  | header :: body ->
    assert(header = dict_header);
    let num_bindings = L.length body in
    let res = Ht.create num_bindings in
    L.iter (fun line ->
        let i, j, mean_sigmas = 
          Scanf.sscanf line "%d-%d:%s" (fun i j rest -> (i, j, rest)) in
        let mean_sigmas = S.split_on_char ',' mean_sigmas in
        let num_dists = L.length mean_sigmas in
        let dists = A.make num_dists { mu = nan; sigma = nan } in
        L.iteri (fun i mean_sigma ->
            dists.(i) <- dist_of_string mean_sigma
          ) mean_sigmas;
        Ht.add res (i, j) dists
      ) body;
    res

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
              -i <filename>: input SMILES file\n  \
              (molecules w/ annotated cut bonds as output by \n  \
              fasmifra_fragment.py)\n  \
              -o <filename>: output file for generated molecules\n  \
              (imcompatible w/ -of)\n  \
              -n <int>: number of molecules to generate\n  \
              [-mu <float>]: average score for all fragments\n  \
              (initial guess; for 1st TS iteration only)\n  \
              [-sigma <float>]: standard deviation for all fragments\n  \
              (initial guess; for all TS iterations)\n  \
              [-of <filename>]: output SMILES fragments to file\n  \
              (not necessarily canonical SMILES; incompatible w/ -o)\n  \
              [-ig <filename>]: load (fragments') gaussians from file\n  \
              [-og <filename>]: output gaussians to file\n  \
              [-pcb]: Preserve Cut Bonds (PCB) in output file\n  \
              [-f]: overwrite existing indexed fragments cache file\n  \
              [-s|--seed <int>]: RNG seed (for repeatable results\n  \
              w/ same input file)\n  \
              [--scores <filename>]: tab-separated name score file\n  \
              (molecule names and order must match the input SMILES file;\n  \
              for Thompson Sampling)\n  \
              [--deep-smiles]: input/output molecules in DeepSMILES\n  \
              no-rings format\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let n = CLI.get_int_def ["-n"] args 0 in
  let input_frags_fn = CLI.get_string ["-i"] args in
  let maybe_scores_fn = CLI.get_string_opt ["--scores"] args in
  let output_fn = CLI.get_string_def ["-o"] args "/dev/null" in
  let maybe_frags_out_fn = CLI.get_string_opt ["-of"] args in
  let maybe_in_gauss_fn = CLI.get_string_opt ["-ig"] args in
  (* don't read/write to same frags dict. file: we want to see
     the evolution of distributions *)
  let maybe_out_gauss_fn = CLI.get_string_opt ["-og"] args in
  let force = CLI.get_set_bool ["-f"] args in
  let use_deep_smiles = CLI.get_set_bool ["--deep-smiles"] args in
  let preserve_cut_bonds = CLI.get_set_bool ["-pcb"] args in
  let maybe_mu = CLI.get_float_opt ["-mu"] args in
  let maybe_sigma = CLI.get_float_opt ["-sigma"] args in
  let get_rng, rng = match CLI.get_int_opt ["-s";"--seed"] args with
    | None -> ((fun x -> x),
               BatRandom.State.make_self_init ())
    | Some seed -> ((fun x -> Random.State.split x),
                    BatRandom.State.make [|seed|]) in
  CLI.finalize (); (* ------------ CLI parsing ----------------------------- *)
  (if Option.is_some maybe_frags_out_fn && output_fn <> "/dev/null" then
     let () = Log.fatal "use either -o (most users) or -of" in
     exit 1
  );
  let use_TS = Option.is_some maybe_mu || Option.is_some maybe_sigma in
  let ij2dists, dists_out_fn =
    match (maybe_in_gauss_fn, maybe_out_gauss_fn) with
    | (None, None) -> (Ht.create 0, "/dev/null")
    | (_, None) | (None, _) ->
      let () = Log.fatal "provide -ig and -og" in
      exit 1
    | (Some in_fn, Some out_fn) ->
      if in_fn = out_fn then
        let () = Log.fatal "-og would overwrite -ig" in
        exit 1
      else
        (load_gaussians in_fn, out_fn) in
  let choose_frag =
    if use_TS then
      thompson_max ij2dists
    else
      uniform_random in
  let assemble =
    if use_deep_smiles then
      if preserve_cut_bonds then
        (Log.fatal "DeepSMILES: unsupported -pcb";
         exit 1)
      else
        assemble_deepsmiles_fragments
    else (* use SMILES *)
    if preserve_cut_bonds then
      assemble_smiles_fragments_PCB
    else
      assemble_smiles_fragments in
  (match maybe_scores_fn with
   | None -> ()
   | Some scores_fn ->
     let () = Log.info "reading scores from %s" scores_fn in
     let smi_name_scores = load_scores input_frags_fn scores_fn in
     match maybe_sigma with
     | None ->
       let () = Log.fatal "--scores requires -sigma" in
       exit 1
     | Some sigma ->
       (* global variance *)
       let s2 = sigma *. sigma in
       let () = Log.info "updating gaussians" in
       update_many_gaussians s2 ij2dists smi_name_scores;
       let () = Log.info "writing updated frags dict. to %s" dists_out_fn in
       (* save_fragments_dict frag2can_smi_id dists dists_out_fn *)
       ()
  );
  Log.info "indexing fragments";
  let seed_fragments, frags_ht =
    let res = load_indexed_fragments maybe_frags_out_fn force input_frags_fn in
    cache_indexed_fragments force input_frags_fn res;
    res in
  Log.info "seed_frags: %d; attach_types: %d"
    (A.length seed_fragments) (Ht.length frags_ht);
  LO.with_out_file output_fn (fun out ->
      let i = ref 0 in
      while !i < n do
        try
          let rng' = get_rng rng in
          let tokens, frag_ids = assemble choose_frag rng' seed_fragments frags_ht in
          fprintf_tokens out tokens;
          (match frag_ids with
           | [] -> fprintf out "\tgenmol_%d\n" !i
           | _ -> fprintf out "\t%s\n" (string_of_frag_ids frag_ids)
          );
          incr i
        with Too_many_rings -> () (* skip it *)
      done
    );
  let stop = Unix.gettimeofday () in
  let dt = stop -. start in
  Log.info "rate: %.2f molecule/s" ((float n) /. dt)

let () = main ()
