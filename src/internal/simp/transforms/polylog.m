(* ====================================================================== *)
(* polylog.m -- FullSimplify identities for the polylogarithm.              *)
(*                                                                          *)
(* Loaded lazily when PolyLog appears in the input. Special values          *)
(* (PolyLog[2,1] -> Pi^2/6, PolyLog[s,0] -> 0, PolyLog[1,z] -> -Log[1-z])   *)
(* are already handled by the kernel; the two-argument dilogarithm          *)
(* functional equations below are the genuine additions.                    *)
(* ====================================================================== *)

(* Euler reflection:  Li2[z] + Li2[1-z] == Pi^2/6 - Log[z] Log[1-z]. *)
RegisterTransforms[PolyLog, {
    Function[e, Replace[e,
        PolyLog[2, z_] + PolyLog[2, w_] /; w === 1 - z :>
            Pi^2/6 - Log[z] Log[1 - z]]],
    (* Landen / duplication:  Li2[z] + Li2[-z] == (1/2) Li2[z^2]. *)
    Function[e, Replace[e,
        PolyLog[2, z_] + PolyLog[2, nz_] /; nz === -z :>
            PolyLog[2, z^2]/2]]
}];
