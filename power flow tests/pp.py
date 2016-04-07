import pkgutil, importlib, inspect
import pypower

# this should probably be turned into modules
# pp.functions and pp.constants, then users could use
# from pp.functions import ___ or pp.functions.___()
# and from pp.constants import *

def import_functions():
    """Import the same-named functions from all pypower.* modules (if the exist).
    
    pypower uses a Matlab-style configuration for files, where each function is 
    defined in its own file. In Matlab, you can access each function just by specifying
    its name, e.g., case30(). But in Python you have to use the awkward equivalent
    from pypower.case30 import case30; case30(). 
    
    import_pypower_functions() adds all these functions to the caller's local namespace,
    equivalent to "from pypower.<function> import <function>".
    
    This is normally called during initialization of the pp module, in which case
    the functions appear as pp.<function>(). But this can also be called from another
    module in which case all the functions will be added to the caller's namespace
    (mostly useful for interactive testing.) 
    
    Alternatively, you can use "from pp import *" or "from pp import case30"."""

    # get a reference to the local variables of the caller (these are the global module
    # variables if this is called from the top level of a module or an interactive session)
    namespace = inspect.stack()[1][0].f_locals    # 1=previous caller; 0=frame object
    try:
        for importer, modname, ispkg in pkgutil.iter_modules(pypower.__path__):
            if not modname.startswith("idx_"):
                # add the same-named function from the module to the caller's local namespace
                # (equivalent to "from pypower.<function> import <function>")
                module = importlib.import_module("pypower." + modname)
                if hasattr(module, modname):
                    namespace[modname] = getattr(module, modname)
        # TODO: patch makePTDF and makeBdc to use bus indexes as bus numbers?
        # i.e., make a copy of the various inputs and renumber the buses 
        # according to their index in the bus list. Or maybe that can be
        # added to these functions in matpower instead?

    finally:
        # avoid creating object cycles
        del namespace
                    


def define_constants(namespace=None):
    """Add all constants defined in the idx_* modules to the specified namespace (e.g.,
    globals()), or to the caller's local namespace if not specified.

    These constants are useful for specifying indexes into the standard arrays for buses,
    branches, etc.
    
    This is equivalent to running "from pypower.idx_<component> import *" for all types 
    of component.
    """
    if namespace is None:
        # get a reference to the local variables of the caller 
        # (which are also the globals in an interactive session or main thread of a module)
        namespace = inspect.stack()[1][0].f_locals    # 1=previous caller; 0=frame object
    try:
        for importer, modname, ispkg in pkgutil.iter_modules(pypower.__path__):
            # import constants from all pypower.idx_* modules
            if modname.startswith("idx_"):
                module = importlib.import_module("pypower." + modname)
                # copy all attributes that don't start with underscores to the appropriate namespace
                for attrib in dir(module):
                    if not attrib.startswith("_"):
                        namespace[attrib] = getattr(module, attrib)
    finally:
        # avoid creating object cycles
        del namespace

import_functions()

