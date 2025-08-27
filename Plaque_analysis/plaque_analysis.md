## Ilastik processing
Ilastik segmentation model and code for processing segmentation results for the paper. Please refer to the [Ilastik documentation](https://www.ilastik.org/) for details on how the program works, and to our paper to learn how it was used there.

### Ilastik model
The "folder contains the compressed classifier file "PlaqueAssay.ilp.zip" and a folder with the training images (10 .tif files). To use the classifier, first download and unzip it in a directory. You'll also need to have the associated training images in a directory. Upon loading the classifier, it will prompt a message that it can't find the images (see the image below). You'll need to manually give the paths to them, and then this new path will be saved for when you re-open the classifier.

<img width="422" height="191" alt="image" src="https://github.com/user-attachments/assets/304c421a-0cb0-47dc-a628-7c75b2779e30" />

## Segmentation data analysis
A set of Jupyter notebooks to process, curate, and quantify the resulting Ilastik segmentations. Runing the first cell of each notebook will install the necessary libraries. Alternatively, you can run `pip install requirements_segmentation.txt` to install them, including Jupyter. We suggest to create first an environment with a manager (Octavio is partial to [miniforge](https://github.com/conda-forge/miniforge)) so there won't be compatibility issues with other software.

### Segmentation processing
The file "watershed.ipynb" contains the code to process the Ilastik segmentations (which only separate the image in background and plaques) into instance segmentations of plaques. We also include example images. Below is an example of the output from running the code and plotting the processed segmentation:

<img width="2100" height="1500" alt="Watershed" src="https://github.com/user-attachments/assets/b8cf0743-8531-4260-bfa3-726f7b051f86" />

### Segmentation curation

The file "plaque_selector.ipynb" contains the code to curate the instance segmentations (obtained with watershed.ipynb). The code will open each image in separate windows where then you can click on **the plaques you wish to keep**. When done, it will create a new image with only the selected plaques.

### Feature extraction 

The file "plaque_quantifier.ipynb" takes the instance segmentation masks (curated or not) and builds a data frame containing the plaque area and intensity quantifications.
