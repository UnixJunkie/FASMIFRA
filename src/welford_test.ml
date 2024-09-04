
(* exercise the Welford module from the CLI
   output mean and stdev of samples read from given file *)

open Printf

module L = BatList
module LO = Line_oriented

let main () =
  if Array.length Sys.argv = 1 then
    (eprintf "usage: %s <FILE>\n\
              FILE must contain one float sample per line\n" Sys.argv.(0);
     exit 1)
  else
    let input_fn = Sys.argv.(1) in
    let samples = LO.map input_fn float_of_string in
    match samples with
    | [] -> assert false
    | x :: xs ->
      let w = L.fold Welford.update (Welford.create x) xs in
      Printf.printf "%.3f\t%.3f\n" (Welford.mean w) (sqrt (Welford.sample_s2 w))

let () = main ()
