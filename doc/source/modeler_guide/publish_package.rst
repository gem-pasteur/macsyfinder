.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2022 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _publish_package:

*************************
Publishing/sharing models
*************************


 .. _writing_model_package:


Writing your own macsy-model package
====================================

The whole package structure and the corresponding files are described in the section :ref:`package_structure`. It requires five different
types of files to be complete:

* a `metadata.yml` file (mandatory)
* a `README.md` file (mandatory)
* a `LICENSE` file (optional but **HIGHLY** recommended)
* a `model_conf.xml` file (optional)
* macsy-models definition(s) within a `definitions` folder (mandatory)
* HMM profiles within a `profiles` folder (mandatory)



Sharing your models
===================

If you want to share your models you can create a :ref:`macsy-model package <model_package>` in your github repository. 
Several steps are needed to publish your model:

1. Check the **validity** of your package with the ``macsydata check`` command.
   You have to run it from within the folder containing your package files. 
   It will report:

   * everything is clear: `macsydata` displays the next step totake to publish the package

   * warning: it means that the package could be improved.
   
   It is better to fix it if you can, but you can also proceed to *Step 2*

   * error: the package is not ready to be published as is. You have to fix the errors before you go to *Step 2*.

2. Create a **tag**, and submit a **pull request** to the https://github.com/macsy-models organization.
   This step is **very important**: without a tag, there is no package.
   `macsydata check` only tagged packages.
   It is also the duty of the model provider to setup a tag with the same name as the version in the `metadata.yml` file.
   It is **Mandatory** to follow a versioning scheme described here:

        * https://www.python.org/dev/peps/pep-0440/#public-version-identifiers
        * https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/specification.html#standard-versioning-schemes

   If your package is in version *2.0.1* the tag must be `2.0.1`.
   The version or tag must **NOT** start with letter as `v2.0.1` or `my_package-2.0.1`.

   .. warning::

        Check that the tag match with the version defined in `metadata.yml`.
        To avoid inadvertent mistake place the script below in `.git/hooks/` directory.
        Check that the hook is well named pre-push and it is executable (`chmod 755 .git/hooks/pre-push`)
        This script check if you push a tag and if the tag match the version in metadata.yml
        If it does not match it prevent the push.

        .. literalinclude:: ../_static/code/pre-push
           :language: shell

        :download:`pre-push <../_static/code/pre-push>` .


3. When your pull request (PR) is accepted, the model package becomes automatically available to the community through the `macsydata` tool.

If you don't want to submit a PR you can provide the tag release tarball (tar.gz) as is to your collaborators.
This archive will also be usable with the `macsydata` tool.

.. note:: 

    ``macsydata check``
    checks the syntax of the package, but it does not publish anything.
    It just warns you if something is wrong with the package.
    Every model provider should check its own package before publishing it.
    The package publication is done by the `git push` and the `pull request`.

Examples of ``macsydata check`` outputs:


Your package is syntactically correct:

.. code-block:: text

    macsydata check tests/data/models/test_model_package/
    Checking 'test_model_package' package structure
    Checking 'test_model_package' metadata_path
    Checking 'test_model_package' Model definitions
    Models Parsing
    Definitions are consistent
    Checking 'test_model_package' model configuration
    There is no model configuration for package test_model_package.
    If everyone were like you, I'd be out of business
    To push the models in organization:
            cd tests/data/models/test_model_package
    Transform the models into a git repository
            git init .
            git add .
            git commit -m 'initial commit'
    add a remote repository to host the models
    for instance if you want to add the models to 'macsy-models'
            git remote add origin https://github.com/macsy-models/
            git tag 1.0b2
            git push --tags


You received some warnings: 

.. code-block:: text

    macsydata check tests/data/models/Model_w_conf/
    Checking 'Model_w_conf' package structure
    Checking 'Model_w_conf' metadata_path
    Checking 'Model_w_conf' Model definitions
    Models Parsing
    Definitions are consistent
    Checking 'Model_w_conf' model configuration
    The package 'Model_w_conf' have not any LICENSE file. May be you have not right to use it.
    The package 'Model_w_conf' have not any README file.
    macsydata says: You're only giving me a partial QA payment?
    I'll take it this time, but I'm not happy.
    I'll be really happy, if you fix warnings above, before to publish these models.

You received some errors:

.. code-block:: text

    macsydata check tests/data/models/TFF-SF/
    Checking 'TFF-SF' package structure
    The package 'TFF-SF' have no 'metadata.yml'.
    Please fix issues above, before publishing these models.
    ValueError
