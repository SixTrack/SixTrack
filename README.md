This branch holds a simple proposal for how to split the sixtrack.s into a giant amount of files (sorted by decks), and build using cmake rather than astuce.

## Instructions

 - Check out branch
 - In folder SixTrack, run ``python sixconv.py``
 - Build using cmake, with SixTrack_noastuce as the top level directory
 
To summarize:
```
cd SixTrack
python sixconv.py
mkdir ../build
cd ../build
cmake ../SixTrack_noastuce
make
```
   
