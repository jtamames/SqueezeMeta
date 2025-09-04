#!/usr/bin/env python

from os import environ
from subprocess import run

import yaml
import re
import jinja2

from argparse import ArgumentParser


def main(args):

    excluded_pkgs = {'squeezemeta',
                     'squeezemeta-dev',
                     '_libgcc_mutex',
                     '_openmp_mutex',
                     '_r-mutex',
                     'alsa-lib',
                     'binutils_impl_linux-64',
                     'binutils_linux-64',
                     'bzip2',
                     'c-ares',
                     'ca-certificates',
                     'certifi',
                     'curl',
                     'dbus',
                     'gcc_impl_linux-64',
                     'gcc_linux-64',
                     'gfortran_impl_linux-64',
                     'gfortran_linux-64',
                     'glib',
                     'glib-tools',
                     'icu',
                     'krb5',
                     'lerc',
                     'libcurl',
                     'libdeflate',
                     'libffi',
                     'libpq',
                     'libuuid',
                     'libxcb',
                     'libxkbcommon',
                     'libxml2',
                     'openssl',
                     'python_abi',
                     'sysroot_linux-64',
                     'xz',
                     'zstd',
                     'motulizer',
                     'superpang',
                     'gtdbtk',
                     'speedict',
                     'igraph',
                     '_x86_64-microarch-level',
                     'kernel-headers_linux-64',
                     'ld_impl_linux-64'
                    }

    ### Read the template and separate the jinja definitions for the yaml content
    template_str = open(args.template)
    jinja_headers = {}
    template_yaml = ''
    jinja_val_buff = ''
    for line in template_str:
        if line.startswith('{%'):
            var, val = line.strip().split('=', 1)
            var = var.strip().split('set ')[-1]
            val = val.lstrip()
            if line.endswith('%}\n'): # single line definition
                jinja_headers[var] = val.rstrip(' %}\n')
            else:
                jinja_val_buff += f'{val}\n' # store in a buffer to allow for multi-line jinja2 var definitions
        elif jinja_val_buff:
            if line.strip().endswith('%}'): # multi line definition, ending in this line
                jinja_val_buff += line.rstrip(' %}\n')
                jinja_headers[var] = jinja_val_buff
                jinja_val_buff = ''
            else: # multi line definition, continued
                jinja_val_buff += line
        else:
            template_yaml += line.replace('{{', '\'{{').replace('}}','}}\'') # quote jinja variables to allow parsing yaml without interpolating them
 
    ### Fix the jinja headers if needed
    jinja_headers['name'] = '"squeezemeta"' if not args.dev else '"squeezemeta-dev"'
    if args.git_tag:
        jinja_headers['git_tag'] = f'"{args.git_tag}"'
    else:
        jinja_headers['git_tag'] = '"master"' if not args.dev else '"develop"'

    ### Parse the yaml section of the template
    template = yaml.safe_load(template_yaml)

    ### Replace package versions and/or add new packages from the current conda environment
    extra_pkgs = {}
    o = run([environ['CONDA_EXE'], 'list'], capture_output=True).stdout.decode()
    for line in o.split('\n'):
        sort_keys=Falseline = line.strip()
        if not line or line.startswith('#'):
            continue
        line = [f for f in line.split(' ') if f]
        if len(line)==3:
            pkg, version, build = line
            channel = ''
        elif len(line)==4:
            pkg, version, build, channel = line
        extra_pkgs[pkg] = {'version': version, 'build': build, 'channel': channel}

    run_reqs = template['outputs'][0]['requirements']['run']
    corrected_run_reqs = []
    added_pkgs = set()
    # Correct versions of the packages in the template to the versions in the current environment, if required
    #  Ignore packages with '{{' in their version, as those will be taken care of by jinja
    for pkg_str in run_reqs:
        pkg = re.split('>|=| ', pkg_str)[0]
        if not args.keep_template_versions and pkg in extra_pkgs and '{{' not in pkg_str:
            newver = extra_pkgs[pkg]['version']
            pkg_str = f'{pkg}=={newver}'
        corrected_run_reqs.append( pkg_str )
        added_pkgs.add(pkg)
    # Add new packages if required
    if not args.ignore_extra_dependencies:
        for pkg, info in extra_pkgs.items():
            version = info['version']
            if pkg in added_pkgs or pkg in excluded_pkgs or info['channel'] == 'pypi':
                continue
            pkg_str = f'{pkg}=={version}'
            corrected_run_reqs.append( pkg_str )

    ### Fix YAML
    if args.env:
        template = {'name': template['package']['name'] + '_' + template['package']['version'],
                    'channels': ['conda-forge', 'bioconda', 'nodefaults'],
                    'dependencies': corrected_run_reqs}
    else:
        template['outputs'][0]['requirements']['run'] = corrected_run_reqs
        del template['about'] # since the yaml parser mangles it, we remove it from here and add it manually later
 
    ### Print output
    out = ''

    # Add the jinja definitions
    for var, val in jinja_headers.items():
        if args.env and var == 'git_tag':
            continue
        out += '{% ' + f'set {var} = {val}' + ' %}\n'
    out += '\n'

    # Add the yaml content

    class MyDumper(yaml.Dumper):
        # We need to create a custom dumper to get the indentation as conda-build likes
        # https://stackoverflow.com/questions/25108581/python-yaml-dump-bad-indentation
        def increase_indent(self, flow=False, indentless=False):
            return super(MyDumper, self).increase_indent(flow, False)

    out_yaml = yaml.dump(template, Dumper=MyDumper, sort_keys=False, default_flow_style=False)
    out_yaml = out_yaml.replace('\'{{', '{{').replace('}}\'', '}}')
    out += out_yaml + '\n'

    # Add the about section manually if needed 
    if not args.env:
        reached_about = False
        for line in open(args.template): # Add the "about" section manually
            if line.startswith('about:'):
                reached_about = True
            if reached_about:
                out += line

    # Render jinja2 template if working with an env file, since conda does not support jinja2 in those
    if args.env:
        out = jinja2.Template(out).render().strip()
    
    # Finally print it to stdout
    print(out)


def parse_args():
    parser = ArgumentParser(description='Create a conda meta.yaml file from a template specifying the package versions used in the current active environment')
    parser.add_argument('--template', help = 'Template meta.yaml file for the dev package', type = str, default = 'template.meta.yaml')
    parser.add_argument('--git-tag', type = str, help = 'Git tag from which to pull SqueezeMeta. Defaults to `master` if `--dev` is not specified, otherwise `develop`')
    parser.add_argument('--dev', action = 'store_true', help = 'Target version is a dev version')
    parser.add_argument('--ignore-extra-dependencies', action = 'store_true', help = 'Do not include dependencies not originally specified in the template')
    parser.add_argument('--keep-template-versions', action = 'store_true', help = 'Keep the package versions originally specified in the template')
    parser.add_argument('--env', action = 'store_true', help = 'Create environment file instead of build file')
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())
