use cxx_build::CFG;
use std::path::Path;
use std::env;
use std::process::Command;

fn main() {
    // build and link slepc_wrapper 
    let boost_root =  env::var("BOOST_ROOT").expect("Environment Variable 'BOOST_ROOT' is not defined!");
    CFG.exported_header_dirs.push(Path::new(&boost_root));
    cxx_build::bridge("./src/slepc_wrapper.rs")
        .file("./cpp_src/slepc_wrapper.cc")
        .flag_if_supported("-std=c++17")
        .flag_if_supported("-lrt")
        .compile("eigensolver");

    // build eigensolver (using makefile)
    Command::new("make").current_dir("./cpp_src/gep_module/").spawn().expect("Error running Makefile for GEP solver!");

    // mark reruns
    println!("cargo:rerun-if-changed=/src/slepc_wrapper.rs");
    println!("cargo:rerun-if-changed=/cpp_src/slepc_wrapper.cc");
    println!("cargo:rerun-if-changed=/cpp_src/slepc_wrapper.h");
    println!("cargo:rerun-if-changed=/cpp_src/gep_module/solve_gep.cpp");
}