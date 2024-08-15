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

(* list all fragments in a SMILES w/ annotated cut bonds *)
let list_fragments smi =
  let seed_frags = ref [] in
  let frags_ht = Ht.create 11 in
  let rewritten = rewrite_paren_cut_bond_smiles smi in
  let tokens = tokenize_full rewritten in
  fragment seed_frags frags_ht tokens;
  Ht.fold (fun _i_j frags acc ->
      L.rev_append frags acc
    ) frags_ht !seed_frags

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
       )
  );
  (seeds, ht)

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

(* a Gaussian distribution *)
type dist = { mu: float;
              sigma: float }

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

(* default fragment sampling policy for training-set distribution matching *)
let uniform_random_choice rng a =
  A.unsafe_get a (BatRandom.State.int rng (A.length a))

(* maximization by Thompson sampling policy *)
let thompson_sample_max rng all_dists frag2_can_smi_id frags =
  (* retrieve the distribution corresponding to each fragment *)
  let dist_frags =
    A.map (fun frag ->
        let _can_smi, frag_id = Ht.find frag2_can_smi_id frag in
        (all_dists.(frag_id), frag)
      ) frags in
  (* play slot machines in Las Vegas *)
  let sample_frags =
    A.map (fun (dist, frag) ->
        (gauss rng dist, frag)
      ) dist_frags in
  let max_score = A.max (A.map fst sample_frags) in
  (* if several have max_score, choose one at random *)
  let candidates = A.filter (fun (score, _frag) -> score = max_score) sample_frags in
  let _max_score, frag = uniform_random_choice rng candidates in
  frag
      
let assemble_smiles_fragments rng seeds branches =
  let frag_count = ref 0 in
  let ht = Ht.create 97 in
  let seed_frag =
    let chosen = uniform_random_choice rng seeds in
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
          let chosen = uniform_random_choice rng possible_branches in
          (* Log.error "chosen branch: %s" (string_of_tokens chosen); *)
          rev_renumber_ring_closures ht frag_count chosen in
        (* Log.error "renumbered branch: %s" (string_of_tokens (L.rev branch)); *)
        (* typed cut bond discarded here *)
        loop acc (L.rev_append branch xs)
      | _  -> loop (x :: acc) xs
  in
  loop [] seed_frag

(* like assemble_smiles_fragments, but "Preserve Cut Bonds" (PCB).
   To output generated molecules w/ cut bonds preserved; so that
   generated molecules do not need to be fragmented later on *)
let assemble_smiles_fragments_PCB rng seeds branches =
  let frag_count = ref 0 in
  let ht = Ht.create 97 in
  let seed_frag =
    let chosen = uniform_random_choice rng seeds in
    L.rev (rev_renumber_ring_closures ht frag_count chosen) in
  let rec loop acc tokens = match tokens with
    | [] -> L.rev acc
    | x :: xs ->
      match x with
      | Cut_bond (i, j) ->
        let possible_branches = Ht.find branches (i, j) in
        let branch =
          let chosen = uniform_random_choice rng possible_branches in
          rev_renumber_ring_closures ht frag_count chosen in
        (* preserve cut bond [x] here *)
        loop (x :: acc) (L.rev_append branch xs)
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
        let branch = uniform_random_choice rng possible_branches in
        (* typed cut bond discarded here *)
        loop acc (L.append branch xs)
      | _  -> loop (x :: acc) xs
  in
  let seed_frag = uniform_random_choice rng seeds in
  loop [] seed_frag

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

let square x =
  x *. x

let stddev (mean: float) (l: float list): float =
  sqrt (L.favg (L.map (fun x -> square (x -. mean)) l))

(* count the number of unique fragment ids in [fn] *)
let num_ids_in_frags_dict fn =
  let lines = LO.lines_of_file fn in
  match lines with
  | [] -> assert false
  | header :: other_lines ->
    (* enforce expected format *)
    assert(header = "smi\tcan_smi\tid\tmean\tstddev");
    let ids =
      L.fold_left (fun acc line ->
          try Scanf.sscanf line "%s@\t%s@\t%d\t%f\t%f"
                (fun _smi _can_smi id _mu _sigma ->
                   ISet.add id acc)
          with exn ->
            (Log.fatal "Fasmifra.load_frag_ids_from_frags_dict: \
                        cannot parse: %s" line;
             raise exn)
        ) ISet.empty other_lines in
    assert(ISet.min_elt ids = 0);
    let n = ISet.cardinal ids in
    assert(ISet.max_elt ids = n - 1);
    Log.info "%s: %d unique fragments" fn n;
    n

(* fragments dictionary header line *)
let dict_header = "#smi\tcan_smi\tid\tmean\tstddev"

(* load in a file created by fasmifra_frag_dict.py *)
let load_fragments_dict maybe_init_dist fn =
  let n = num_ids_in_frags_dict fn in
  let frag2can_smi_id = Ht.create n in
  let use_global_dist, init_dist = match maybe_init_dist with
    | None ->
      (* each line is supposed to have non NaN mu and sigma then *)
      (false, { mu = nan; sigma = nan})
    | Some x ->
      (assert (not (Float.is_nan x.mu) &&
               not (Float.is_nan x.sigma));
       (true, x)) in
  let dists = A.make n { mu = nan; sigma = nan} in
  LO.with_in_file fn (fun input ->
      let header = input_line input in
      (* enforce expected format *)
      assert(header = dict_header);
      while true do
        let line = input_line input in
        try Scanf.sscanf line "%s@\t%s@\t%d\t%f\t%f"
              (fun smi can_smi id mu sigma ->
                 let frag =
                   tokenize_full (rewrite_paren_cut_bond_smiles smi) in
                 (* the canonical SMILES is unused at run-time
                  * (but useful to check the fragments dictionary)
                  * we keep it so that we can output a complete dictionary *)
                 Ht.add frag2can_smi_id frag (can_smi, id);
                 if Float.is_nan mu || Float.is_nan sigma then
                   (assert use_global_dist;
                    dists.(id) <- init_dist)
                 else
                   dists.(id) <- { mu; sigma }
              )
        with exn ->
          (Log.fatal "Fasmifra.load_fragments_dict: cannot parse: %s" line;
           raise exn)
      done
    );
  (frag2can_smi_id, dists)

(* save fragments dictionary w/ updated distributions to file *)
let save_fragments_dict frag2can_smi_id dists fn =
  (* sort fragments by (id, smi) *)
  let n = Ht.length frag2can_smi_id in
  let arr = A.make n ("", "", -1) in
  let i = ref 0 in
  Ht.iter (fun frag (can_smi, id) ->
      let frag_smi = string_of_tokens frag in
      arr.(!i) <- (frag_smi, can_smi, id);
      incr i
    ) frag2can_smi_id;
  A.sort (fun (smi0, _can_smi0, id0) (smi1, _can_smi1, id1) ->
      if id0 < id1 then
        -1
      else if id0 > id1 then
        1
      else (* id0 = id1 *)
        BatString.compare smi0 smi1
    ) arr;
  LO.with_out_file fn (fun out ->
      fprintf out "%s\n" dict_header;
      A.iter (fun (smi, can_smi, id) ->
          let dist = dists.(id) in
          fprintf out "%s\t%s\t%d\t%f\t%f\n" smi can_smi id dist.mu dist.sigma
        ) arr
    )

(* update the Gaussian score distribution for each fragment
   [s2] is the initial guess for the variance of all fragments *)
let update_gaussians s2 frag2can_smi_id dists smi_name_scores: unit =
  L.iter (fun (smi, _mol_name, score) ->
      L.iter (fun frag ->
          (* all fragments are supposed to be in the frags dictionary *)
          let _can_smi, id = Ht.find frag2can_smi_id frag in
          dists.(id) <- update_gaussian dists.(id) s2 score
        ) (list_fragments smi)
    ) smi_name_scores

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
              -n <int>: number of molecules to generate\n  \
              [-mu <float>]: average score for all fragments\n  \
              (initial guess)\n  \
              [-sigma <float>]: standard deviation for all fragments\n  \
              (initial guess)\n  \
              [-of <filename>]: output fragments to SMILES file\n  \
              (not necessarily canonical ones)\n  \
              [-ifd <filename>]: load fragments dictionary file\n  \
              (first run must also use -mu and -sigma)\n  \
              [-ofd <filename>]: output fragments dictionary\n  \
              [-pcb]: Preserve Cut Bonds (PCB) in output file\n  \
              [-f]: overwrite existing indexed fragments cache file\n  \
              [-s|--seed <int>]: RNG seed (for repeatable results\n  \
              w/ same input file)\n  \
              [--scores <filename>]: tab-separated name score file\n  \
              (molecule names and order must match the input SMILES file;\n  \
              for Thompson sampling)\n  \
              [--deep-smiles]: input/output molecules in DeepSMILES\n  \
              no-rings format\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let n = CLI.get_int ["-n"] args in
  let input_frags_fn = CLI.get_string ["-i"] args in
  let maybe_scores_fn = CLI.get_string_opt ["--scores"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let maybe_frags_out_fn = CLI.get_string_opt ["-of"] args in
  let maybe_frags_dict_in_fn = CLI.get_string_opt ["-ifd"] args in
  let maybe_frags_dict_out_fn = CLI.get_string_opt ["-ofd"] args in
  let force = CLI.get_set_bool ["-f"] args in
  let use_deep_smiles = CLI.get_set_bool ["--deep-smiles"] args in
  let preserve_cut_bonds = CLI.get_set_bool ["-pcb"] args in
  let maybe_mu = CLI.get_float_opt ["-mu"] args in
  let maybe_sigma = CLI.get_float_opt ["-sigma"] args in
  let assemble =
    if use_deep_smiles then
      if preserve_cut_bonds then
        (Log.fatal "PCB mode for DeepSMILES unsupported yet";
         exit 1)
      else
        assemble_deepsmiles_fragments
    else (* use SMILES *)
    if preserve_cut_bonds then
      assemble_smiles_fragments_PCB
    else
      assemble_smiles_fragments in
  let get_rng, rng = match CLI.get_int_opt ["-s";"--seed"] args with
    | None -> ((fun x -> x),
               BatRandom.State.make_self_init ())
    | Some seed -> ((fun x -> Random.State.split x),
                    BatRandom.State.make [|seed|]) in
  CLI.finalize (); (* ------------ CLI parsing ---------------- *)
  let maybe_init_dist = match (maybe_mu, maybe_sigma) with
    | (Some mu, Some sigma) -> Some { mu; sigma }
    | _ -> None in
  let frag2can_smi_id, dists, frags_dict_out_fn =
    match (maybe_frags_dict_in_fn, maybe_frags_dict_out_fn) with
    | (None, None) -> (Ht.create 0, [||], "/dev/null")
    | (Some _, None) ->
      let () = Log.fatal "-ifd without -ofd" in
      exit 1
    | (None, Some _) ->
      let () = Log.fatal "-ofd without -ifd" in
      exit 1
    | (Some in_fn, Some out_fn) ->
      if in_fn = out_fn then
        let () = Log.fatal "-ofd would overwrite -ifd" in
        exit 1
      else
        let ht, dists = load_fragments_dict maybe_init_dist in_fn in
        (ht, dists, out_fn) in
  let smi_name_scores = match maybe_scores_fn with
    | None -> []
    | Some scores_fn -> load_scores input_frags_fn scores_fn in
  (match smi_name_scores with
   | [] -> ()
   | _ ->
     match maybe_init_dist with
     | None ->
       let () = Log.fatal "--scores requires -mu and -sigma" in
       exit 1
     | Some init_dist ->
       (* global variance *)
       let s2 = init_dist.sigma *. init_dist.sigma in
       let () = Log.info "updating gaussians" in
       update_gaussians s2 frag2can_smi_id dists smi_name_scores;
       let () = Log.info "writing updated frags dict. to %s" frags_dict_out_fn in
       save_fragments_dict frag2can_smi_id dists frags_dict_out_fn
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
          let tokens = assemble (get_rng rng) seed_fragments frags_ht in
          fprintf_tokens out tokens;
          fprintf out "\tgenmol_%d\n" !i;
          incr i
        with Too_many_rings -> () (* skip it *)
      done
    );
  let stop = Unix.gettimeofday () in
  let dt = stop -. start in
  Log.info "rate: %.2f molecule/s" ((float n) /. dt)

let () = main ()
