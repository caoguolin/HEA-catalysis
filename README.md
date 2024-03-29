# **HEA-catalysis**
The code and data for a "HEA-catalysis" project

The code contains the "High throughput calculations framework" and "Co-adsorption model"

## **High throughput calculations framework**
The high-throughput calculations framework is primarily composed of two parts: batch modeling of HEA surfaces and DFT high-throughput calculations. For the batch modeling of HEA surfaces, we firstly constructed a Pt(111) surface as model system. Then, we randomly replaced Pt atoms with the nine noble metal elements based on  special quasi-random structures (SQS) method, and generated 100 HEA surfaces. Based on these surfaces, we finally constructed a total of 8170 adsorption systems. For the DFT high-throughput calculations, we built a task dispatcher for automatic management of the calculations, containing check submission, file uploading & downloading, job submission & monitoring three subroutines:

_**Check submission**_. First of all, the dispatcher will check whether the calculation tasks have been submitted or not. If not, the dispatcher will upload files and construct the submission scripts from scratch. This will from a queue composed of tasks to be executed. Otherwise, the dispatcher will recover the existing queue and collect results from those tasks that have finished.

_**File uploading & downloading**_. The uploading files include four VASP input files (INCAR, POSCAR, POTCAR and KPOINTS). The downloading files are the VASP calculated energy files and optimized structure files of the adsorption systems. Files uploading and downloading have been performed with os and shutil modules in Python. Finally, the calculated results were transmitted and stored in JavaScript Object Notation format to the local device.

_**Job submission & monitoring**_. After uploading needed files, the dispatcher would generate scripts onto sever for executing DFT calculation tasks with a homemade job scheduling system. Then, the job monitoring subroutine monitors the “state” of the tasks during the whole calculation processes and would correct the errors and resubmit it on-the-fly.

## **Co-adsorption model**
Firstly, a 277 Å×277 Å HEA surface was constructed based on SQS method. Then, the adsorption energy was predicted for each active site (on-top site for O adsorption fcc hollow site for OH adsorption) on the HEA surface with our proposed descriptor, and subsequently sort all the active sites according to the values of adsorption energy. The O/OH adsorbs step-by-step on the HEA surface according the following rule: i) the starting adsorption site would be the one with the lowest adsorption energy; ii) for the O adsorbed fcc hollow site, its three neighboring on-top site should avoid the OH adsorption, preventing the so-called blocking effect; iii) the adsorption would continue until no more free surface sites remained.
