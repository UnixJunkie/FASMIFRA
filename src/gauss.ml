
module Log = Dolog.Log

(* a Gaussian distribution *)
type t = { mu: float; (* mean *)
           s2: float (* sigma^2 (variance) *) }

let create mu s2 =
  { mu; s2 }

let to_string d =
  (* think about m+/-s, but shorter *)
  Printf.sprintf "%g/%g" d.mu d.s2

let of_string s =
  Scanf.sscanf s "%f/%f" (fun mu s2 ->
      if s2 = 0.0 then
        (Log.fatal "Gauss.of_string: s2=0"; exit 1)
      else
        { mu; s2 }
    )

let almost_one = Float.pred 1.0
(* enforce that the highest float we can draw is <1.0 *)
let _ = assert(almost_one < 1.0)
let pi = 4.0 *. (atan 1.0)
let two_pi = 2.0 *. pi

(* [sample mu sigma] get one float from the normal distribution
   with mean=mu and stddev=sigma
   a = cos(2*pi*x) * sqrt(-2*log(1-y))
   b = sin(2*pi*x) * sqrt(-2*log(1-y)) (b is ignored below)
   cf. Python's documentation of random.gauss function *)
let sample rng dist =
  dist.mu +. ((sqrt dist.s2) *.
              (cos (two_pi *. (Random.State.float rng almost_one)) *.
               sqrt (-2.0 *. log (1.0 -. (Random.State.float rng almost_one)))))

(* sample [n] times *)
let samples rng dist n =
  Array.init n (fun _i -> sample rng dist)
