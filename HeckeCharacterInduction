// Code by Luca Dall'Ava and Aleksander Horawa. January 2024.
//
// Magma code for associating a pair of Hecke character to the starting CM wt 1 modular form g is a wt 1 modular form with CM by K
//chi:=Character(g);


//Return vector of norms of ideals in "decomposition". Here decomposition:=Factorisation(cond_norm*OK).

Norms_Vector:=function(decomposition);
	
	n:=#decomposition;

	V:=[];

	for i in [1..n] do
		V := Append(V,Norm(decomposition[i][1]));
	end for;

	return V;

end function;

//Find indexes of the vectors of norms. NOT OPTIMAL! I'm computing all combinations and then throwing away the bad ones.

IndexesProduct:=function(vector_norms,target);

list:=[* *];

n:=#vector_norms;

for i in [1..n] do

	V := Subsets( {1 .. n}, i);

	for v in V do
		y:=1;
		for x in v do
			y:=y*vector_norms[x];
	    end for;
	    if y eq target then
	    	list cat:=[* IndexedSetToSequence(SetToIndexedSet(v)) *];
	    end if;
	end for;
end for;
return list;
end function;

// Routine for finding the associated character.
//Sloppy bound! FIX!

FindPsiCouple:=function(g,I,OK : grp:=Elements(HeckeCharacterGroup(I)), N:=0, bound:=Level(g) )
 	
 	G:=grp;
 	index:=0;
	for psi in G do
		index:=index+1;
		if (#G eq 2) then
			return G;
			break;
		end if;
		n:=N;
		while true do
			n:=n+1;
			p:=NthPrime(n);
			if p ge bound then
				break;
			end if;
			P:=p*OK;
			if IsInert(p, OK) then
				if Coefficient(g,p) ne 0 then
					return false;
				end if;
			elif IsSplit(p,OK) then
				fact:=Factorisation(P);
				psi_vals:=0;
				for pp in fact do
					psi_vals := psi_vals + psi(pp[1]);
				end for;
				if Coefficient(g,p) ne psi_vals then
					G:=Remove(G,index);
					index:=index-1; //remove and then we need to shift back
					break;
				end if;
			end if;
		end while;
	end for;
	if #G eq 2 then
		return G;
	else   // If we are not done, we increase the bound ad start with the given G + given N.
		return $$(g,I,OK : grp:=G, N:=n-1, bound:=p^2); 
	end if;
end function;

// Main
// 

HeckeCharacterInduction:=function(g,K);

N:=Level(g);


d:=Discriminant(K);

cond_norm:=N/d;
cond_norm:=Integers()!cond_norm;
divisors:=Factorisation(cond_norm);

OK:=MaximalOrder(K);
fact:=Factorisation(cond_norm*OK);

FACT:=Norms_Vector(fact);

FACT_indexes:=IndexesProduct(FACT,Abs(cond_norm));
	
JJ:=[* *];

for index in FACT_indexes do
	I:=1*OK;
	for i in [1..#index] do
		J:=fact[index[i]][1];
		I:=I*J;
	end for;

	JJ:=Append(JJ, [* FindPsiCouple(g,I,OK), I *]);
	
end for;

return JJ;

end function;

// Last routine to compute the tensor product L-function and output 
//the analytic ranks computed up to the d-th derivative


AnalyticRank_E_K_psi_g:=function(E,K,g:Derivative:=0,Sign)

CHARS:=HeckeCharacterInduction(g,K);
RES:=[* *];
for i in [1..#CHARS] do
	
	psi_g:=CHARS[i][1][1]/CHARS[i][1][2];
	L_chi:=LSeries(Ch);

	LEKch:=TensorProduct(LSeries(E,K), L_chi : Sign:=Sign);

	res:=[];
	
	for j in [0..Derivative] do
		res:=Append(res, Evaluate(LEKch,1:Derivative:=j));
	end for;
	RES:=Append(RES,[*psi_g, res*]);
end for;
return RES;

end function;


//// EXAMPLE
/*

R<x>:=PolynomialRing(Integers());
pol:=x^2-x+2;
K:=NumberField(pol);
M:=ModularForms(Gamma1(175),1);
MD:=DihedralForms(M);
for V in MD do 		
	if Conductor(V[1]) eq 7 then
		C:=V;               
	end if;
end for;
g:=C[2][1];

E_eq:=x^3-947*x+11214;
E:=EllipticCurve(E_eq);

HeckeCharacterInduction(g,K);

AnalyticRank_E_K_psi_g(E,K,g:Derivative:=1, Sign:=1);

*/
