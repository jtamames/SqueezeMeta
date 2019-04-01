"""
Part of the SqueezeMeta distribution. 25/03/2018 Original version, (c) Fernando Puente-Sánchez, CNB-CSIC.
python utilities for working with SqueezeMeta results
"""


def parse_conf_file(project_path):
    perlVars = {}
    for line in open('{}/SqueezeMeta_conf.pl'.format(project_path)):
        line = line.rsplit('#',1)[0] # Remove comment strings.
        if line.startswith('$'): # Is this a var definition?
            var, value = [x.strip(' \'\"') for x in line.strip().strip(';').split('=',1)]
            perlVars[var] = value

    ### Define this bc it's funny to parse perl code with python.
    def perl_string_interpolation(string):
        if '$' in string:
            for var in perlVars:
                if var in string and not '\\{}'.format(var) in string: # The var is in the string, and the $ is not escaped like "\$"
                    string = string.replace(var, perl_string_interpolation(perlVars[var])) # Recursive interpolation.
        return string


    ### Back to work. Interpolate all the strings.
    for var, value in perlVars.items():
        perlVars[var] = perl_string_interpolation(value)

    return perlVars
