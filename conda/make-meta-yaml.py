#!/usr/bin/env python

from os import environ
from subprocess import run

from argparse import ArgumentParser


def main(args):

    pkgs = {}
    excludedPkgs = {'squeezemeta',
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
                   }

    name = 'squeezemeta' if not args.dev else 'squeezemeta-dev'
    if args.git_tag:
        git_tag = args.git_tag
    else:
        git_tag = 'master' if not args.dev else 'develop'

    if not args.keep_template_versions:
        o = run([environ['CONDA_EXE'], 'list'], capture_output=True).stdout.decode()
        for line in o.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            line = [f for f in line.split(' ') if f]
            if len(line)==3:
                pkg, version, build = line
                channel = ''
            elif len(line)==4:
                pkg, version, build, channel = line
            pkgs[pkg] = {'version': version, 'build': build, 'channel': channel}
    


    with open(args.template) as infile: #input meta.yaml
        section = ''
        quotes1 = 0
        quotes2 = 0
        addedRunNotice = False
        addedExtraPkgs = False
        addedPkgs = set()
        lead_white = '  ' if args.env else '    '
        for line in infile:
            line = line.rstrip() # leave leading whitespace
            quotes1 += line.count('"')
            sline = line.strip()
            if not sline.startswith('-') and ':' in sline and not quotes1%2 and not quotes2%2: # this will fail if multi-line strings are created using ' or contain inner "
                section = sline.split(':')[0]
                #print(section)
            if section == 'name' and sline: # this won't work if name or version are defined in several lines
                if not args.env:
                    line = f'  name: {name}'
                else:
                    line = f'name: {name}_{args.version}'
                    line += '\nchannels:'
                    line += '\n  - conda-forge'
                    line += '\n  - bioconda'
                    line += '\n  - anaconda'
            if section == 'version' and sline:
                line = f'  version: "{args.version}"'
            if section == 'git_tag' and sline:
                line = f'  git_tag: {git_tag}'
            if section == 'description' and 'squeezemeta-dev' in line:
                line = line.replace('squeezemeta-dev', name)
            if sline == 'run:' and args.env:
                line = 'dependencies:'
            if section == 'run' and args.env and 'pin_compatible' in line:
                line =  line.replace('{{ ','').replace('{{','').replace(' }}','').replace('}} ','')
                line =  line.replace('pin_compatible', '').replace('( ', '').replace('(','').replace(' )','').replace(')', '')
                line =  line.replace('\'','').replace('"','')
                sline = sline.replace('{{ ','').replace('{{','').replace(' }}','').replace('}} ','')
                sline = sline.replace('pin_compatible', '').replace('( ', '').replace('(','').replace(' )','').replace(')', '')
                sline = sline.replace('\'','').replace('"','')

            if not args.keep_template_versions:
                if section == 'run':
                    if not sline:
                        continue
                    if sline.startswith('-'):
                        if not addedRunNotice:
                            print(f'{lead_white}# Packages originally specified in the template')
                            addedRunNotice = True
                        pkg = sline.replace('>','=').split('=')[0].replace('- ','')
                        if 'numpy' in pkg:
                            pkg = 'numpy'
                        if pkg in pkgs:
                            if pkg not in excludedPkgs:
                                version = pkgs[pkg]['version']
                                line = f'{lead_white}- {pkg}=={version}'
                            else:
                                line = f'{lead_white}- {pkg}'
                            addedPkgs.add(pkg)
                        else:
                            line = f'{lead_white}{line.strip()}'
                        
                if not args.ignore_extra_dependencies and section != 'run' and addedRunNotice and not addedExtraPkgs:
                    print(f'{lead_white}# Extra dependencies')
                    for pkg in pkgs:
                        if pkg not in addedPkgs and pkg not in excludedPkgs:
                            version = pkgs[pkg]['version']
                            print(f'{lead_white}- {pkg}=={version}')
                    print()
                    addedExtraPkgs = True


            quotes2 += line.count('"')

            if args.env and section not in {'name', 'run'}:
                continue

            print(line)


def parse_args():
    parser = ArgumentParser(description='Create a conda meta.yaml file from a template specifying the package versions used in the current active environment')
    parser.add_argument('version', help = 'Target package version')
    parser.add_argument('--template', help = 'Template meta.yaml file for the dev package', type = str, default = 'template.meta.yaml')
    parser.add_argument('--git-tag', type = str, help = 'Git tag from which to pull SqueezeMeta. Defaults to `master` if `--dev` is not specified, otherwise `develop`')
    parser.add_argument('--dev', action = 'store_true', help = 'Target version is a dev version')
    parser.add_argument('--ignore-extra-dependencies', action = 'store_true', help = 'Do not include dependencies not originally specified in the template')
    parser.add_argument('--keep-template-versions', action = 'store_true', help = 'Keep the package versions originally specified in the template')
    parser.add_argument('--env', action = 'store_true', help = 'Create environment file instead of build file')
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())
