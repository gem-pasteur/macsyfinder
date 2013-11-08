


import os


#########################
# Utilities
#########################
def which(name, flags=os.X_OK):
    """Search PATH for executable files with the given name."""
    result = []
    path = os.environ.get('PATH', None)
    if path is None:
        return []
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
            break
    return result
