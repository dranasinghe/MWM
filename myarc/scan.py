#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import datetime
from arkane.statmech import determine_qm_software
import numpy as np
import matplotlib.pyplot as plt

from myarc import settings
from myarc.main import ARC
from myarc.job.ssh import SSHClient, write_file
from myarc.common import read_yaml_file, save_yaml_file
from myarc.parser import parse_xyz_from_file
from myarc.job.job import Job
from myarc.settings import servers, check_status_command, submit_command, submit_filename, delete_command

def disply_scan(project_directory, job_name):
    species_label = find_species_by_job_name(project_directory, job_name)
    log = determine_qm_software(os.path.join(project_directory,'calcs','Species',species_label,job_name,'output.out'))
    print (os.path.join(project_directory,'calcs','Species',species_label,job_name,'output.out'))
    v_list, angle = log.loadScanEnergies()
    v_list = np.array(v_list, np.float64)
    v_list *= 0.001
    plt.plot(angle, v_list,'.-')
    plt.show()

def scan_jobs_in_directory(project_directory):
    species_dict = {}
    species_root_dir = os.path.join(project_directory, "calcs", "Species")
    for species in os.listdir(species_root_dir):
        species_dir = os.path.join(species_root_dir, species)
        if os.path.isdir(species_dir):
            species_dict[species] = [job_type for job_type in os.listdir(species_dir) \
                                     if
                                     (os.path.isdir(os.path.join(species_dir, job_type)) and job_type != "conformers")]
    return species_dict


def find_species_by_job_name(project_directory, job_name):
    species_dict = scan_jobs_in_directory(project_directory)
    # Find the corresponding species
    for species_label, jobs in species_dict.iteritems():
        if job_name in jobs:
            return species_label
    else:
        raise Exception("Specified job not found")


def find_job_path(JOB, project_director, species_label, job_name, servers, server):
    local_path = os.path.join(project_director, "calcs", "Species", species_label, job_name)
    species_name_for_remote_path = species_label.replace('(', '_').replace(')', '_')
    remote_path = os.path.join('/home', servers[server]["un"], 'runs', 'ARC_Projects', JOB.project,
                               species_name_for_remote_path, job_name)
    return local_path, remote_path


def find_scan_info_from_input_file(local_path, input_filename):
    scan_trsh = "\n"
    pivot = None
    local_input_path = os.path.join(local_path, input_filename)
    with open(local_input_path, "r") as f:
        for line in f.readlines():
            if re.search("D[\s\d]+[\s\d]+[\s\d]+[\s\d]+[\s]+S", line.strip()):
                pivot = [int(line.split(" ")[2]), int(line.split(" ")[3])]
                scan = [int(line.split(" ")[1]), int(line.split(" ")[2]), int(line.split(" ")[3]),
                        int(line.split(" ")[4])]
                scan_res = float(line.split(" ")[-1])
            if re.search("D[\s\d]+[\s\d]+[\s\d]+[\s\d]+[\s]+F", line.strip()):
                scan_trsh += re.search("D[\s\d]+[\s\d]+[\s\d]+[\s\d]+[\s]+F", line.strip()).group() + '\n'
    if scan_trsh == "\n":
        scan_trsh = ''
    if pivot:
        return (pivot, scan, scan_res, scan_trsh)
    else:
        return (None, None, None, None)


def get_job_memory(JOB, servers, server):
    memory = JOB.memory
    cpus = servers[server].get('cpus', 8)  # set to 8 by default
    if software == 'molpro':
        memory = memory * 128 / cpus
    if software == 'terachem':
        memory = memory * 128 / cpus
    elif software == 'gaussian':
        memory = memory * 1000
    elif software == 'orca':
        memory = memory * 1000 / cpus
    elif software == 'qchem':
        memory = memory  # dummy
    elif software == 'gromacs':
        memory = memory  # dummy
    return memory


def generate_running_scan_job_dict(input_dict):
    running_jobs = {}
    for species, jobs in input_dict["running_jobs"].iteritems():
        running_jobs[species] = {}
        for job in jobs:
            if job["job_type"] == "scan":
                running_jobs[species][job["job_name"]] = job["pivots"]
        if not running_jobs[species]:
            running_jobs.pop(species, None)
    return running_jobs


def read_job_num(output_path):
    job_name = output_path.split("/")[-2]
    job_num = int(job_name.split('_')[-1][1:])
    return job_num


def compare_pivot(pivot1, pivot2):
    if pivot1 == pivot2 or pivot1 == [pivot2[-1], pivot2[0]]:
        return True
    else:
        return False


def obtain_scan_job_status(input_dict):
    rotor_status = {}
    running_scan_jobs = generate_running_scan_job_dict(input_dict)

    for species in input_dict["species"]:
        rotor_status[species["label"]] = {}
        for rotor in species["rotors_dict"].itervalues():
            # 1. Check whether running:
            if species["label"] in running_scan_jobs.keys():
                for job_name, pivot in running_scan_jobs[species["label"]].iteritems():
                    if compare_pivot(pivot, rotor["pivots"]):
                        rotor_status[species["label"]][str(rotor["pivots"])] = ("Running", job_name)
                        break
            if str(rotor["pivots"]) in rotor_status[species["label"]].keys():
                continue

            # 2. Check if it is a job converged before any corrections
            job_num = 0
            if rotor["scan_path"]:
                job_num = read_job_num(rotor["scan_path"])
                if "opt" in input_dict["output"][species["label"]].keys():
                    job_opt_num = read_job_num(input_dict["output"][species["label"]]["opt"])
                elif "composite" in input_dict["output"][species["label"]].keys():
                    job_opt_num = read_job_num(input_dict["output"][species["label"]]["composite"])
                else:
                    rotor_status[species["label"]][str(rotor["pivots"])] = (
                    "Error: unknown reason 1", "scan_a" + str(job_num))
                    continue
                # Scan can not be earlier than opt
                if job_num < job_opt_num:
                    pivot_from_scan = rotor["scan"][1:3]
                    if compare_pivot(pivot_from_scan, rotor["pivots"]):
                        rotor_status[species["label"]][str(rotor["pivots"])] = (
                        "Error: Converged before comformer correction", "scan_a" + str(job_num))
                        continue
                    else:
                        rotor_status[species["label"]][str(rotor["pivots"])] = (
                        "Error: scan changed to {0} given pivot {1}".format(rotor["scan"], rotor["pivots"]),
                        "scan_a" + str(job_num))
                        continue
                # Check symmetry
                if "symmetry" in rotor.keys():
                    if rotor["success"] and rotor["symmetry"]:
                        rotor_status[species["label"]][str(rotor["pivots"])] = (
                        "Successfully converged", "scan_a" + str(job_num))
                        continue
                    else:
                        rotor_status[species["label"]][str(rotor["pivots"])] = (
                        "Error: unknown reason 1", "scan_a" + str(job_num))
                        continue
                elif rotor["success"]:
                    rotor_status[species["label"]][str(rotor["pivots"])] = (
                    "Error: End job without troubleshooting", "scan_a" + str(job_num))
                    continue

                # job converged, but the range is too large
                if "larger" in rotor["invalidation_reason"]:
                    rotor_status[species["label"]][str(rotor["pivots"])] = (
                    "Converged but invalid", "scan_a" + str(job_num))
                    continue
                # job not converged
                if "crash" in rotor["invalidation_reason"]:
                    # job not converged and no output in local
                    if not os.path.isfile(rotor["scan_path"]):
                        rotor_status[species["label"]][str(rotor["pivots"])] = (
                        "Error: Scan never start or unexpected termination", "scan_a" + str(job_num))
                        continue
                    with open(rotor["scan_path"]) as f:
                        lines = f.readlines()
                    if not lines:
                        rotor_status[species["label"]][str(rotor["pivots"])] = (
                        "Error: Empty output caused by unexpected termination", "scan_a" + str(job_num))
                        continue
                    for line in lines[-1:-20:-1]:
                        if 'Normal termination of Gaussian.py' in line:
                            rotor_status[species["label"]][str(rotor["pivots"])] = (
                            "Error: unknown reason 2", "scan_a" + str(job_num))
                            break
                    else:
                        rotor_status[species["label"]][str(rotor["pivots"])] = (
                        "Error: ESS not converged", "scan_a" + str(job_num))
                        continue
                    continue
                else:
                    rotor_status[species["label"]][str(rotor["pivots"])] = (
                    "Error: unknown reason 3", "scan_a" + str(job_num))
                    continue
            else:
                if species["label"] in input_dict["running_jobs"].keys():
                    for job in input_dict["running_jobs"][species["label"]]:
                        if job["job_type"] == "composite" or job["job_type"] == "opt":
                            rotor_status[species["label"]][str(rotor["pivots"])] = ("Running preliminary job", "N/A")
                            break
                elif "opt" not in input_dict["output"].keys() or "composite" not in input_dict["output"].keys():
                    rotor_status[species["label"]][str(rotor["pivots"])] = (
                    "Running preliminary job or conformer jobs fail", "N/A")
                else:
                    rotor_status[species["label"]][str(rotor["pivots"])] = (
                    "Error: unknown reason 4", "scan_a" + str(job_num))
    return rotor_status