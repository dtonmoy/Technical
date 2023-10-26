### Zip and Split files/folders into multiple segments


**Step 1:** Go to the terminal and navigate to the directory where you will generate the zipped files:
```
cd /directory/of/interest/
```
**Step 2:** Next, do the following:
```
zip -r -s 1900m myzip.zip myfolder
```
**Step 3:** The documentation:
1. "myzip.zip" is the name of the zipped file that you want to create
2. "myfolder" is the name of the folder you want to zip and split
3. "1900m" is the size of the zipped files. You can change this according to your preference