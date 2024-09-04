
(* exercise the Gauss module from the CLI *)

open Printf

let main () =
  if Array.length Sys.argv = 1 then
    (eprintf "usage: %s mu stdev n\n\
              <mu:float> <stdev:float> <n:int>\n" Sys.argv.(0);
     exit 1)
  else
    let mu    = float_of_string Sys.argv.(1) in
    let sigma = float_of_string Sys.argv.(2) in
    let n     = int_of_string   Sys.argv.(3) in
    let rng = BatRandom.State.make [|314159265|] in
    let dist = Gauss.create mu (sigma *. sigma) in
    for _i = 1 to n do
      Printf.printf "%f\n" (Gauss.sample rng dist)
    done

let () = main ()
