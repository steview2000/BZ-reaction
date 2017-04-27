# BZ-reaction 

This repositiory contains the C-code for simulating the BZ reaction. The different programs are:

## BZ-single 
Just a single oregonator

## BZ-couple 
Two coupled oregonators

## BZ-couple 
Three coupled oregonators


## Installation

### For Mac OS X
Necessary dependencies:
git (probably already instealled)
gcc (probably already installed)

The whole code can be installed as follows:

Installing libtiff with homebrew:
	
	`brew install libtiff`

Clone this repository and enter BZ-reaction directory:
	
	`git clone https://github.com/steview2000/BZ-reaction`
	`cd BZ-reaction`


Now compile first the files important for writing .tiff files:
	
	`cd analysisPic`
	`make`

Now compile BZ-couple
	
	```bash
	cd ../BZ-couple
	make
	```
