# Plaque_analysis
Ilastik segmentation model and code for processing segmentation results for the paper. Please refer to the [Ilastik documentation](https://www.ilastik.org/) for details on how the program works, and to our paper to learn how it was used there.

## Ilastik model
The "folder contains the compressed classifier file "PlaqueAssay.ilp.zip" and a folder with the training images (10 .tif files). To use the classifier, first download and unzip it in a directory. You'll also need to have the associated training images in a directory. Upon loading the classifier, it will prompt a message that it can't find the images (see the image below). You'll need to manually give the paths to them, and then this new path will be saved for when you re-open the classifier.

<img width="422" height="191" alt="image" src="https://github.com/user-attachments/assets/304c421a-0cb0-47dc-a628-7c75b2779e30" />

## Segmentation processing

The file ".ipynb" contains the code to process the Ilastik segmentations (which only separate the image in background and plaques) into instance segmentations of plaques. Below is an example of the output from running the code and plotting the processed segmentation:

<img width="2700" height="1500" alt="Watershed_example" src="https://github.com/user-attachments/assets/4bd126bb-6875-47a1-8812-b70c008ff88b" />
