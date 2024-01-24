Given an elliptic curve E, a dihedral weight 1 modular form g, and a quadratic (CM) field K:

1) "HeckeCharacterInduction" outputs couples of Hecke characters over K (with the correspongin conductor ideal) such that the Galois representation of g is isomorphic to the induction from K of any of those characters.

  INPUTS: 
  entries: g ModFormElt, K NmbFld
  parameters: O: OrdNmbFld Default:=MaximalOrder(K);
  

3) "AnalyticRank_E_K_psi" is utilized to compute the analytic rank of L(E/K,psi), for psi varying through the quotients of the characters in the couples obtained with "HeckeCharacterInduction".
