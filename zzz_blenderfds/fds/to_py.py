"""BlenderFDS, tokenize FDS file in a readable notation"""

import re

nl_re = re.compile(r"""
    (?P<namelist>   # namelist, group "namelist"
        ^&                # starting ampersand after newline (re.MULTILINE)
        (?P<label>[a-zA-Z][a-zA-Z0-9_]+?)  # namelist label, not greedy
        [,\s\t]+          # one or more separators of any kind    
        (?P<params>       # namelist params, protect strings, no &
            (?: '[^']*?' | "[^"]*?" | [^'"&] )*?  # zero or more groups, not greedy
        )
        [,\s\t]*          # zero or more separators of any kind    
        /                 # anything outside &.../ is a comment and is ignored
    )
    """, re.VERBOSE | re.MULTILINE | re.DOTALL)

def _extract_namelists(text):
    """Return a list of multiline namelists strings from an fds file"""
    return re.findall(nl_re, text)

param_re = re.compile(r"""
    (?P<label>[a-zA-Z][a-zA-Z0-9_\(\):,]+?)  # parameter label w bounds, not greedy
    [\s\t]*           # zero or more spaces    
    =                 # an equal sign
    [\s\t]*           # zero or more spaces    
    (?P<fds_value>    # the value group, protect strings
        (?: '[^']*?' | "[^"]*?" | [^'"] )+?  # one or more groups, not empty
    )
    (?=               # stop the previous value match when it is followed by
        [,\s\t]+      # one or more separators of any kind
        [a-zA-Z][a-zA-Z0-9_\(\):,]+  # another parameter label (same definition as before)
        [\s\t]*       # zero or more spaces
        =             # an equal sign
        |             # or
        $             # the end of the string
    )
    """, re.VERBOSE | re.DOTALL)  # no MULTILINE, so that $ is the end of the file

def _extract_params(text):
    """Return a list of parameters"""
    return re.findall(param_re, text)

def _eval_param(text):
    """Eval text to the corresponding Py value"""
    # Remove newlines
    text = ' '.join(text.splitlines())
    # Get logicals, or eval
    if text.upper() in ('T','.TRUE.'):
        return True
    elif text.upper() in ('F','.FALSE.'):
        return False
    else:
        return eval(text)

def tokenize(text):
    """Parse and tokenize fds text.
    Input:  "&OBST ID='Hello' XB=1,2,3,4,5,6 /"
    Output: (("&OBST ID='Hello' XB=1,2,3,4,5,6 /", "OBST", {'ID': 'Hello', ...}), ...)
    """
    tokens = list()
    for nl in _extract_namelists(text):
        # nl[0]: original namelist
        # nl[1]: fds label
        params = dict()
        for par in _extract_params(nl[2]):
            # par[0]: fds label
            # par[1]: fds value
            try:
                params[par[0]] = (_eval_param(par[1]), par[1])  # {label: (value, fds_value), }
            except:
                raise Exception(f'Cannot evaluate parameter: {par[0]}={par[1]}')
        # tokens = ("fds_label", {label: (value, fds_value), label: (value, fds_value), ...}, "original namelist"), ...
        tokens.append((nl[1], params, nl[0]))
    return tokens

if __name__ == "__main__":
    import sys
    if not sys.argv: exit()
    print("BFDS fds.to_py.tokenize:", sys.argv[1])
    with open(sys.argv[1], 'r') as f:
        fds_file = f.read()
    print('\n'.join(str(v) for v in tokenize(fds_file)))

