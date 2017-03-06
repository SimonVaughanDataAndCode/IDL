;+
; NAME:
; Index_Clumps
; PURPOSE:
; Isolate groups of consecutive subscripts.
; Procedure returns two arrays: Clump_Beg and Clump_End
; which specify all the consecutive runs in the input array SubScripts.
; The consecutive subscripts are then given by the Loop:
;   for j=0,n_elements(Clump_Beg)-1 do begin
;     SubClump = SubScripts( Clump_Beg(j):Clump_End(j) )
; CALLING:
; Index_Clumps, SubScripts, Clump_Beg, Clump_End, NOSINGLE= , COUNT=
; INPUT:
; SubScripts = array of numbers.
; KEYWORDS:
; COUNT = the number of clumps found (output).
; /NOSINGLE : ignore trivial clumps of one subscript.
; OUTPUTS:
; Clump_Beg = indices of clump beginings.
; Clump_End = indices of clump endings.
; PROCEDURE:
; Apply the where function.
; HISTORY:
; Written: Frank Varosi NASA/GSFC 1991.
; FV 1997, fixed for the case of just one subscript.
;-

Pro Index_Clumps, SubScripts, Clump_Beg, Clump_End, NOSINGLE=nosin, COUNT=Nclump

  Nsub = N_elements( SubScripts )
  Clump_Beg=[0]
  Clump_End=[Nsub-1]  ;default values.

  if (Nsub LE 0) then begin
    print,"syntax:  Index_Clumps, SubScripts, Clump_Beg, Clump_End"
    print,"keywords:    COUNT=Nclump, /NOSINGLE"
    Nclump = 0
    RETURN
  endif

  if (Nsub EQ 1) then begin
    Nclump = 1
    RETURN
  endif

  Bound = WHERE( ABS( SubScripts[1:*]-SubScripts[*] ) GT 1, Nclump )
  Nclump = Nclump+1
  if (Nclump EQ 1) then RETURN

  Clump_Beg = [0, Bound+1]
  Clump_End = [Bound, Nsub-1]

  IF KEYWORD_SET( nosin ) THEN Begin

    w = WHERE( Clump_End GT Clump_Beg, Nclump )

    IF (Nclump GT 0) Then begin
      Clump_Beg = Clump_Beg[w]
      Clump_End = Clump_End[w]
    endif
  endif
end

