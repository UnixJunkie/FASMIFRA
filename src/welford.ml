
(* Welford algorithm for computing mean and variance on the fly
 * cf. https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance *)

type t = { mu: float; (* mean *)
           sm2: float; (* sum of squared deviations to the mean *)
           n: float } (* number of measurements as a float (for perf.) *)

let create mu =
  { mu; sm2 = 0.0; n = 1.0 }

(* take into account new measurement *)
let update prev x =
  let n = prev.n +. 1.0 in
  let mu = prev.mu +. (x -. prev.mu) /. n in
  let sm2 = prev.sm2 +. (x -. prev.mu) *. (x -. mu) in
  { mu; sm2; n }

let mean x =
  x.mu

(* sample variance (biased) *)
let sample_s2 x =
  x.sm2 /. x.n

(* population variance (unbiased) *)
let pop_s2 x =
  x.sm2 /. (x.n -. 1.0)
