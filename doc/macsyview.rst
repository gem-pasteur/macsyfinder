.. _macsyview:

MacSyView: visualizing MacSyFinder's results!
===============================================

MacSyView is a browser application allowing to visualize MacSyFinder's detected systems. 

****************
MacSyView How To
****************

1. Run MacSyFinder to detect your favorite system!
2. Launch MacSyView
3. Select the .json output file in the output directory of the run
4. Choose the system you want to visualize in the list...
5. ...and here it is! 


.. note::
    The MacSyView application runs everything on the user's computer, even if it uses the technologies of Web browsers. No data are sent out of the user's device.

****************************
Graphical output description
****************************

The content of the system view depends on the type of the input dataset. 

- upper panel: an overview of the effectives of detected components is displayed **for all types** of datasets, per type of components in the system definition. It is a direct representation of how the definition was fulfilled during detection.

- middle panel: only for **ordered datasets**, the detected system is shown in its genomic context, including nearby proteins that were not annotated as system's components.

- lower panel: for **all datasets**, a table containing information on detected components (and eventually nearby proteins for ordered datasets) is displayed. It includes sequence information, and in the case of system's components, Hmmer hit information, and function assigned in the system. 


***************
Technology used
***************




