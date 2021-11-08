*=============================================================================*
* Project title: staggered_did
* Created by: Pietro Santoleri (p.santoleri@santannapisa.it)
* Created on: 19/07/2021
* Last modified on: 06/08/2021
* This is a modified version of a code written by Borusyak and Burgherr 
*=============================================================================*

// Clear Memory
clear all


// Set options 
*=============================================================================*

// Download user-written packages (0 for no; 1 for yes)
global downloads 0

// Run DID script (0 for no; 1 for yes)
global run_did 1


// Use included packages
*=============================================================================*

cap adopath - PERSONAL
cap adopath - PLUS
cap adopath - SITE
cap adopath - OLDPLACE
adopath + "stata_packages"
net set ado "stata_packages"

// Download packages 
if $downloads == 1 {
	// Install  Packages
	ssc install blindschemes, replace
	ssc install estout, replace	
	ssc install ftools, replace
	ssc install reghdfe, replace
	ssc install did2s, replace
	* net install did2s, from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/") replace
	ssc install did_imputation, replace
	ssc install event_plot, replace
	ssc install did_multiplegt, replace              
	ssc install eventstudyinteract, replace
	ssc install drdid, replace     // new DID estimator
	ssc install flexpaneldid, replace      // new DID estimator
	ssc install xtevent, replace        // new DID estimator
	ssc install eventdd, replace     // new DID estimator
	ssc install csdid, replace				
	* net install csdid, from ("https://raw.githubusercontent.com/friosavila/csdid_drdid/main/code/") replace
	* 注意：guithub连不上，下面的三个包不能加载；如github可以连上，就应该加载它们
	net install github, from("https://haghish.github.io/github/") replace
	* github install lsun20/eventstudyinteract, replace
	github install joshbleiberg/stackedev, replace
	
	* 使用staggered_stata包需要先安装recall包
	github install haghish/rcall, version("2.5.0")
	github install jonathandroth/staggered_stata, replace       // new DID estimator
}

// Set your file paths, These are relative to where the project file is saved. 
*=============================================================================*

global data_path "data/data_for_analysis"
global raw_data_path "data/raw_data" 
global temp_path "temp" 

global script_path "scripts" 
global output "output" 
global log_path "scripts/logs" 

global restricted_data_raw "restricted_raw_data" 
global restricted_data_analysis "restricted_data_for_analysis"

// Version of stata
version 17

// Specify Screen Width for log files
set linesize 255

// Set font type
graph set window fontface "Roboto"

// Set Graph Scheme
graph set window fontface default
set scheme plotplainblind

// Allow the screen to move without having to click more
set more off

// Drop everything in mata
matrix drop _all

// Run do files
*=============================================================================*

// Run script
if $run_did == 1 {
	do $script_path/1.staggered_did_analysis.do
}




