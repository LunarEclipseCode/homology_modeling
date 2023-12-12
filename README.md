## Fixes

I saw Victor tried with two new proteins. The issue with AF_AFA0A023IWE3F1 is it's too small sequence, we need something larger at the moment. Also, when choosing a protein, please choose something with high global pLDDT score. (So, when you go to the site, you want to sort by global plDDT score (best to worst). The higher the pLDDT, the more confident AlphaFold was when it predicted the 3D model. The idea behind choosing higher pLDDT score model is if AlphaFold cannot predict it well, our simple model is not doing it. AF_AFA0A023IWE3F1 has very low pLDDT score (about 61) so we are not using it.

For the MA_MAT3VR3956 model, it's a multi chain model, and comparative modeling isn't good enough to do multi-chain stuff.

## Updates

As the BLAST thing takes a while, I went through the setup and now you can run it locally. You can just git pull and all the necessary files are inside the blast-2.15... that folder.

Also, I have looked through bunch of proteins, you can try with this one AF_AFB9K865F1, the sequence alignment for this has bunch of gaps, so it's a bit difficult one.

### More Update and Stuff To Do

I have added a new blast_new file. So, I instead of trying to make blast.py more general by adding more if condition and cases, I was trying to do the multiple sequence alignment thing we mentioned like in 1 slide. In this approach, I am trying to get rid of the handling gaps problem. 

But it's not really multiple sequence alignment, I will explain it in a second. So, in blast.py, in the server search, we used the protein that matches the best with the input sequence. But in blast_new.py, we are using an array of proteins that matchest closest with the input sequence. In the array, the first pdb id matches the best, the second one matches less than the first id and so on.

I am currently working on the following step:

Now, we do sequence alignment of input sequence with all the protein sequences. Now, we have about 200 protein sequences, so I am hoping for example, index 1-50, exactly matches with some protein, then index 40 - 75 matches with some other protein, and so on.... So, we will combine pdb files from different proteins and form the new pdb. Now, even though I am trying to avoid gaps and mismatch as much as possible in this method, there can still be gaps and mismatch but again it should be much less than the blast.py approach. But again, if we are combining multiple pdb files, then the post processing gets more complicated. Also, AlphaFold uses multiple sequence alignment, so if we can get blast_new to work, then we will have a much robust program.


