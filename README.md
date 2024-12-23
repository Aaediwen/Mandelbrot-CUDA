# Mandelbrot-CUDA
A mandelbrot browser using CUDA and Xlib

This was written entirely as an experiment with CUDA programming

build: 

    nvcc mandelbrot.cu -ltiff -lX11 -lm -o mandelbrot

usage:

    ./mandelbrot [Winx=1024] [WinY=768 (max 1080)] 

    q -- Exit
    c -- zoom in
    d -- zoom out
    r -- cuda single xwin
    f -- cuda single tiff
    y -- test xwin
    h -- test tiff
    t -- cpu xwin
    g -- cpu tiff
    ' ' -- preview render
    right click -- cuda xwin
    middle click -- cuda tiff
    left click -- preview render
    scroll -- zoom
