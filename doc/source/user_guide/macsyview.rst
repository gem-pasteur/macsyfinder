.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014  Institut Pasteur, Paris.                           
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _macsyview:

MacSyView: visualizing MacSyFinder's results!
===============================================

MacSyView is a standalone web-browser application to visualize MacSyFinder's detected systems.
MacSyView relies on JSON files outputted by MacSyFinder to display the list of detected systems,
and a detailed view of each system. It allows visualizing the content of systems, their genomic context,
and generates SVG files that can be exported for drawing purpose.

****************
MacSyView How To
****************

1. Run MacSyFinder to detect your favorite system!
2. Launch MacSyView: 

 * Either, **run the wrapper `macsyview`** installed with MacSyFinder's binaries (*i.e.*, `macsyfinder` - for Linux). 
 * Or **open with your web-browser the html page**: /usr/share/macsyview/index.html or /share/macsyview/index.html
   (or in the path specified during installation for data associated with MacSyFinder).

3. Select the .json output file in the output directory of the run
4. Choose the system you want to visualize in the list...
5. ...and here it is! 


.. note::
    The MacSyView application runs everything on the user's computer, even if it uses the technologies of Web browsers.
    No data are sent out of the user's device.

****************************
Graphical output description
****************************

The content of the system view depends on the type of the input dataset. 

 * upper panel: an overview of the effectives of detected components is displayed **for all types** of datasets,
   per type of components in the system definition.
   It is a direct representation of how the definition was fulfilled during detection.
 * middle panel: only for **ordered datasets**, the detected system is shown in its genomic context,
   including nearby proteins that were not annotated as system's components.
 * lower panel: for **all datasets**, a table containing information on detected components
   (and eventually nearby proteins for ordered datasets) is displayed. It includes sequence information,
   and in the case of system's components, Hmmer hit information, and function assigned in the system.


***************
Technology used
***************

MacSyView was coded in Javascript and uses third-party libraries that are all accredited in the COPYRIGHT file
distributed with the MacSyFinder/MacSyView package.

It includes among others:
 
 * the `Raphael library <http://raphaeljs.com/>`_ for systems drawing, 
 * the `Bootstrap library <http://getbootstrap.com/>`_ for HTML design and 
 * the `Mustache library <http://github.com/janl/mustache.js>`_ for HTML templating in Javascript. 
 
The `JQuery <http://jquery.com/>`_, `JQuery-mousewheel <https://github.com/brandonaaron/jquery-mousewheel>`_ and
`Raphael.Export <http://github.com/ElbertF/Raphael.Export>`_ libraries were also used.
It was tested on Chromium and Firefox for Linux, and on Chrome, Firefox and Safari for Mac OS X. 

.. _screenshot:

**********
Screenshot
**********

Here is a view of one of the three systems detected with the example dataset :ref:`presented here <datatest>`:

    .. image:: ../_static/fig_capture.*
     :height: 600px
     :width: 700px 
     :align: left

