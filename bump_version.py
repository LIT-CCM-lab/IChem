#!/usr/bin/python3

import sys
import os
import re
import subprocess
from datetime import datetime

WK_DIRPATH = os.path.dirname(os.path.realpath(__file__))

def get_version_from_cmake(cmake_filepath):
    line_regex = re.compile("^project\(IChem VERSION \d+\.\d+\.\d+\)$")
    version_regex = re.compile("\d+\.\d+\.\d+")

    major, minor, path = None, None, None
    with open(cmake_filepath) as f:
        for line in f:
            re_result = line_regex.match(line)
            if re_result:
                if major is not None:
                    raise ValueError("Several regex matches! Only one is expected.")
                major, minor, path = [int(x) for x in version_regex.search(line).group().split('.')]
                break

    return major, minor, path


def get_version_from_switch(switch_filepath):
    line_regex = re.compile('^#define ICHEM_VERSION "\d+\.\d+\.\d+"$')
    version_regex = re.compile("\d+\.\d+\.\d+")

    major, minor, path = None, None, None
    with open(switch_filepath) as f:
        for line in f:
            re_result = line_regex.match(line)
            if re_result:
                if major is not None:
                    raise ValueError("Several regex matches! Only one is expected.")
                major, minor, path = [int(x) for x in version_regex.search(line).group().split('.')]
                break

    return major, minor, path


def get_release_date(switch_filepath):
    line_regex = re.compile('^#define ICHEM_RELEASE "\d{4}-\d{2}-\d{2}"$')
    release_date_regex = re.compile("\d{4}-\d{2}-\d{2}")

    year, month, day = None, None, None
    with open(switch_filepath) as f:
        for line in f:
            re_result = line_regex.match(line)
            if re_result:
                if year is not None:
                    raise ValueError("Several regex matches! Only one is expected.")

                year, month, day = release_date_regex.search(line).group().split('-')

    return year, month, day

def main(args):
    cmake_filepath = os.path.join(WK_DIRPATH, "CMakeLists.txt")
    if not os.path.isfile(cmake_filepath):
        print("CMakeLists.txt not found!", file=sys.stderr)
        return 1

    switch_filepath = os.path.join(WK_DIRPATH, "source", "headers", "ICTools", "switch.h")
    if not os.path.isfile(switch_filepath):
        print("switch.h not found!", file=sys.stderr)
        return 1

    major, minor, path = get_version_from_cmake(cmake_filepath)
    if major is None:
        print("Major version not found!", file=sys.stderr)
        return 1

    if minor is None:
        print("Minor version not found!", file=sys.stderr)
        return 1

    if path is None:
        print("Path version not found!", file=sys.stderr)
        return 1

    if (major, minor, path) != get_version_from_switch(switch_filepath):
        print("Versions are different!")
        return 1

    year, month, day = get_release_date(switch_filepath)
    if year is None:
        print("Year not found!", file=sys.stderr)
        return 1

    if month is None:
        print("Month not found!", file=sys.stderr)
        return 1

    if day is None:
        print("Day not found!", file=sys.stderr)
        return 1

    if args.increment == "major":
        new_major = major + 1
        new_minor = 0
        new_path = 0
    elif args.increment == "minor":
        new_major = major
        new_minor = minor + 1
        new_path = 0
    elif args.increment == "path":
        new_major = major
        new_minor = minor
        new_path = path + 1
    else:
        print("Wrong increment argument!", file=sys.stderr)

    print(f"Version will be bumped to {new_major}.{new_minor}.{new_path}")
    answer = input("Continue? (yes/no) >")
    while answer not in ('yes', 'no'):
        print("Wrong answer!")
        answer = input("Continue? (yes/no) >")

    if answer == 'no':
        print("Nothing to do. Exiting...")
        return 0

    cmd = [
        "sed", "-i",
        rf"s/VERSION {major}.{minor}.{path}/VERSION {new_major}.{new_minor}.{new_path}/g",
        cmake_filepath
    ]

    process = subprocess.run(cmd, stderr=subprocess.PIPE)

    if len(process.stderr) != 0:
        print("Something went wrong while changing version number in cmake file!", file=sys.stderr)
        return 1

    current_date = datetime.now()
    cmd = [
        "sed", "-i", "-e",
        rf's/ICHEM_VERSION \"{major}.{minor}.{path}\"/ICHEM_VERSION \"{new_major}.{new_minor}.{new_path}\"/g',
        "-e",
        rf's/ICHEM_RELEASE \"{year}-{month}-{day}\"/ICHEM_RELEASE \"{current_date.strftime("%Y-%m-%d")}\"/g',
        switch_filepath
    ]

    process = subprocess.run(cmd, stderr=subprocess.PIPE)

    if len(process.stderr) != 0:
        print("Something went wrong while changing version number in switch file!", file=sys.stderr)
        return 1

    print("Version successfully bumped!")
    print("Please compile the source code in release mode and copy the new version in the 'distrib' directory.")
    print(f"Use this message for git commit: 'Release {new_major}.{new_minor}.{new_path}'.")
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser('Change automatically version number in source code')

    parser.add_argument('-i', '--increment', choices=["major", "minor", "path"], required=True)
    parser.set_defaults(func=main)

    args = parser.parse_args()

    # Run
    status = args.func(args)
    sys.exit(status)


