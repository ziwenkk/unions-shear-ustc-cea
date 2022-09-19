'crossdr72' is the executable file, to run the code, you can enter the following command on a Linux system: './crossdr72 /path/file_name &'

The galaxy data (file_name) should be in a three-column format: ra, dec, redshift.

If you want to change the cosmological parameters or the separation bins, you can modify those parameters in 'crossdr72.f' (line 154 to 171).

After modification, compile the code by running the command: "sh spt.sh".

