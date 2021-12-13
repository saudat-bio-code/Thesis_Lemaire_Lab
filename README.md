This reposityory is organised in the following way:

- Folder dependecy - This folder includes vital packages neccessery for analysis using ciona/crobu genomes
- Folder scripts - subdirectories for analysis in prep negative sets, svm training/testing, ancestry reconstruction 
- Folder data - includes every intermidiated files, such as peak files, svm results, radndom negative sets, etc
- Folder figures - self explainotory, this folder is for plots
- Folder ancestry - this folder is for ancestor genome reconstructions


GitHub changed access options recently from password to toyken/ deplyoed key options.
I configured this repo to be accessible by the same rsa-key that I generated for belem2 login access.
Usually this line of code works well:
     sudo -u salishayeva ssh -T git@github.com -i public_key
     Enter passphrase for key 'public_key': polki
public_key.rsa and public_key files are provided in the dependency/access folder for convenience.

After rsa-key configured it is important to load the repo using ssh not hhttp connection, this command works well

git clone git@github.com:saudat-bio-code/Thesis_Lemaire_Lab.git
#should work without password


![Alt text](/figures/intro.png)
