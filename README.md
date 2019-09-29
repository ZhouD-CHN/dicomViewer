Dicom-RT File Viewer Prototype
===============================
This is a prototype of a Dicom-RT File Viewer. The functions used in the GUI are described in the dicomanal.py. To see the details of the functions, please read the dicomanal.py and the comments in it. This introduciton is focused on the GUI, especially the function of the buttons. Please contact me if you have any suggestions or questions.

****
|Author|Zhou Dejun|
|---|---
|E-mail|dejunzhou@outlook.com
****

![][dicomViewer/ScreenShot.png]

# 1.How to start.
This program is based on Python3, the required Modules are listed in the file requirements.txt. To run this program, put the viewergui.py and dicomanal.py in the same folder, cd to this folder and run the command "pythonw viewergui.py". Please note that the "pythonw" is used here instead of "python" here to avoid the "This program needs access to the screen" problem.

# 2.Open Patient Button
Click this button and select the folder that stores the '.dcm' files of the patient. Then the path of the folder will be displayed on the text box next to it and the image of the first dcm file (sorted by the name of the files) will be displayed.

# 3.Previous Slice Button
Click this button to show the previous image of the dcm file on the file list. If the first slice is being displayed, click this button and the last slice will be displayed.

# 4.Next Slice Button
Click this button to show the next image of the dcm file on the file list. If the last slice is being displayed, click this button and the first slice will be displayed.

# 5.Go Button
Input the number of slice in the textbox before it, then click this button and the required slice will be displayed.

# 6.Show Structure Check Box
If this check box is checked, the sturcture region will be displayed on the image. The legend of the structure will be showed on the left on the image.

# 7.Show Isodose Check Box
If this check box is checked, the dose distribution will be displayed on the image. The legend of the isodose will be showed on the right on the image. Please note that if there is no dose distribution on the image, nothing will be changed if this check box is checked.
