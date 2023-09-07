
### R version control for mac users

You need to follow the following steps for installing and using multiple versions of R in mac:

**Step 1:** Go to the terminal and navigate to the following folder:
```
cd /Library/Frameworks/R.framework/Versions
```

**Step 2:** You will find different version of R installed here (possibly). Now the goal is to make your mac forget about the installed version of R. From the R documentation if you run the following two lines, you will get the list of files that mac need to forget.

```
pkgutil --files org.r-project.x86_64.tcltk
pkgutil --files org.r-project.x86_64.texinfo
```

**Step 3:** Now run the forget command in the directory.
```
pkgutil --forget org.r-project.x86_64.tcltk
pkgutil --forget org.r-project.x86_64.texinfo
```

**Step 4:** Once the pkgutil trickery is done, you can go ahead and download the *.pkg of the R version that you want to install and install it like you normally do. You will find the installed new version of the R in the same directory you previously navigated.

**Step 5:** Install Rswitch which is a nice looking GUI. You are successfully done. Just open Rswitch and change to the version you want to use.