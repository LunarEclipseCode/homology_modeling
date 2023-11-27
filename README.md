## Everything about the code so far.

#### Changes to previous code

* I have changed the code for the part where given input sequence, we use BLAST to find the protein that has the most similar sequence. The previous code was stuck for some reason.
  
  The new code works but don't try with new proteins too often because then the server doesn't respond to too frequent requests. You can use a VPN to bypass this problem though.
  
  Also,  fior this reason, I have added a mapping.txt file that when I am asking the server for a match for the first time, the response will be saved to the mapping file, so that we don't have to ask the server anymore if we need it later.

* Also changed the sequence alignment implementation because the pairwise2 function gives deprecation warning. This took forever :(
  
  
  

#### Newly added stuff

So, now that we have the alignment and we can download the pdb of the template protein, we have to create the 3D model / the PDB file of the input sequence.

First I tried with [AF_AFA0A023GPI8F1](https://www.rcsb.org/structure/AF_AFA0A023GPI8F1) and it happens to that it exactly matches with [4K1Y](https://www.rcsb.org/structure/4K1Y)) because it happens to be that AF_AFA0A023GPI8F1 is one of the 8 identical chains in 4K1Y. The extract_chain function gets one chain from 4K1Y. Now, the RCSB site contains a 3D model generated of AF_AFA0A023GPI8F1 using AlphaFold program and we want to comparehow similar our model is to theirs. However, the RCSB site's model is in CIF format and so cif_to_pdb.py is to convert that CIF file to PDB so that we can compare.

Now, to compare, we use RMS (root mean square) and use only carbon atoms because the carbon atoms form the backbone of a protein. We get 0.6 angstom RMS which means our model is very similar to AlphaFold. Now, you might think, well, if it's exact match in sequence, why don't just extracting one chain give RMS of 0 angstrom. That's because when you get one chain from the system, your are ignoring the inter-chain interaction and so the atom's positions in the edge of the chain will not be accurate in our model.

#### About to be added (I have the idea, currently working on the code)

The next one s a bit more complicated. [AF_AFB5EZH0F1](https://www.rcsb.org/structure/AF_AFB5EZH0F1) almost matches with one chain in #D46. But the key is 'almost' which means there is gaps in the sequence alignment. Now, if there is gaps, how do we predict the atom position when there is a mismatch?

I initially thought about finding a general system but that's too complicated. So, I will code for this specific case and as we later test with more porteins, add more stuff to the code.

So, if you see the sequence alignment in this case, there are multiple ones with the same score but one is clearly better than the rest (the 2nd alignment in sequence.txt file). However, changing scoring function is not that easy in this case. So, I am currently writing a function that filters between the alignments. Basically the idea is if two alignments are exactly same after index = i where i should be relatively small, for now, let's say below 10. Then, for characters before index i, I will favor the alignment that has starts with a gap for input sequence but no gap in template sequence.

Now, before editing a PDB file, we need to understand what different columns mean in a pdb file. Please read the first 2 pages of pdbformat.pdf file in this repo before reading the rest of this file. It took so long to come up with the following idea because there is literally no open source code for this, no paper actually details in depth how they create the 3D model. So, even though it's a re-implementaion project, from here on, it's  orginal work.

Now if we are trying to create the PDB of 2nd alignment in sequence.txt, the one that looks like:

----MTLPKIKHVRAWFIGGATAEKGAGGGDY

MENIMTLPKIKHVRAWFIGGATAEKGAGGGDY

Then, looking at the PDB file of temp.pdb (one chain of 3D46), we can see, we should start from line 34 because the previous lines correspond to M (MET), E (GLU), N(ASN), I(ILM). Now, if we are starting from line 34. So, line 34 will be line 1 in the pdb file of AF_AFB5EZH0F1. However, that means we need to change atom serial number and the residue sequence number. **IF YOU WANT TO START CODING, PLEASE TRY TO DO THIS PART.** Because there is a bit more biology involved part that I have to think about.

Well, if you follow the sequence alignment, you will after the initial gap, the only mismatch is here.

ENG-QT

EN-RQT
So, for G (GLY), it's gap in template, but for next one, template has R (ARG), while input has gap, Now II was thinking if they both have same number of atoms, then we can replace R's atomic co-ordinates with G's and just make other changes. However, that is not the case here. If you look at file aa_format_pdb.txt, you will see R needs 11 atoms in PDB while G is 9 atoms. So, I have to either figure out the peptide bond length between N and G and G and Q or I have to look for a smilar pattern in a different protein and figure out the co-ordinates.

Also, the complicated part is if G's coordinates is changed, that means you cannot directly copy paste template's coordinates for the rest of the sequence, you have to take changes into account and adjust accordingly.


