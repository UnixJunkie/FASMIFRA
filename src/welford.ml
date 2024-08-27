
(* Welford's algorithm for computing mean and variance on the fly
 * cf. https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance *)

module A = BatArray
module Log = Dolog.Log

type t = { mu: float; (* mean *)
           sm2: float; (* sum of squared deviations to the mean *)
           n: float } (* number of measurements as a float (for perf.) *)

(* constructor from a single sample *)
let create mu =
  { mu; sm2 = 0.0; n = 1.0 }

(* special constructor (if you know all internal values) *)
let create_full mu sm2 n =
  { mu; sm2; n }

let to_string w =
  Printf.sprintf "%g/%g:%g" w.mu w.sm2 w.n

let of_string s =
  Scanf.sscanf s "%f/%f:%f" (fun mu sm2 n ->
      if sm2 = 0.0 || n < 2.0 then
        (Log.fatal "Welford.of_string: forbidden params: %s" s; exit 1)
      else
        { mu; sm2; n }
    )

(* update rule (take into account new measurement) *)
let update prev x =
  let n = prev.n +. 1.0 in
  let mu = prev.mu +. (x -. prev.mu) /. n in
  let sm2 = prev.sm2 +. (x -. prev.mu) *. (x -. mu) in
  { mu; sm2; n }

(* constructor from an array of samples *)
let of_samples arr =
  assert(A.length arr >= 2);
  A.fold_lefti (fun acc i x ->
      if i > 0 then
        update acc x
      else (* i = 0 *)
        create x
    ) (create_full 0. 0. 0.) arr

let mean x =
  x.mu

(* sample variance (biased) *)
let sample_s2 x =
  x.sm2 /. x.n

(* population variance (unbiased) *)
let pop_s2 x =
  x.sm2 /. (x.n -. 1.0)
