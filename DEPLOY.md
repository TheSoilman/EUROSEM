# Hosting the EUROSEM web app (free, public)

Target: Render's free web service, built from `Dockerfile`. Not tested on Render or with
Docker yet (the author's sandbox had no Docker daemon or gfortran); do step 1 first.

## What public mode does

`server.py --public` (the Dockerfile sets `EUROSEM_PUBLIC=1`):

* sensitivity analysis is off (routes return 404, the page removes the section)
* each visitor: at most 20 model runs per 10 minutes; everybody together: 40 per minute
  (`--rate-runs`, `--rate-window`, `--global-runs`); over the limit gives a plain message
* listens on 0.0.0.0 and on `$PORT`; `GET /healthz` returns `ok`
* the visitor is the first `X-Forwarded-For` address; a determined user could vary it to dodge
  the per-visitor limit, which is why the all-visitors limit exists
* nothing is stored: each run uses a temporary folder that is deleted afterwards

## 1. Check the container locally (on your Mac, needs Docker Desktop)

The build context is the folder that holds `headless/` and `eurosem_web/`
(`~/Downloads/eurosem`), with `Dockerfile`, `.dockerignore` and `render.yaml` copied into it.

    cd ~/Downloads/eurosem
    docker build -t eurosem .
    docker run --rm -p 8080:8080 eurosem
    # open http://localhost:8080 ; run the example; it should give 26.84 mm and 7878 kg

The build compiles the model with the same flags as `build_exe.sh`. No Docker? Skip to step 2:
Render builds it for you, and the build log shows any compile error.

## 2. GitHub repository

Create a repository (private works on Render) containing exactly:

    Dockerfile  .dockerignore  render.yaml
    headless/stage4_rval/*.for          (the model source; the other stages are not needed)
    eurosem_web/                        (the app; tests and reference may stay, they are ignored)

Decide here whether the model source should be public: a public repository publishes the commented
Fortran. Add a licence file if it is. (Model and guide authorship is yours to state.)

## 3. Render

1. render.com, sign up with GitHub (no card needed for the free plan).
2. New + > Blueprint, choose the repository (or New + > Web Service > Docker, Free plan,
   health check path `/healthz`).
3. Wait for the build. The service gets an address like `https://eurosem-xxxx.onrender.com`.

Free-plan behaviour (render.com/docs/free): sleeps after 15 minutes idle and takes about a minute
to wake; 512 MB RAM and under one CPU; 750 instance hours a month; suspended for the month if
bandwidth runs out (no card on file); "not for production".

## 4. WordPress (personal plan)

Add a page or menu link to the Render address. Say that the first visit after idle time can take
a minute to load. (A personal-plan site generally cannot embed custom apps, so link rather than
embed.)

## Updating

Push to the repository; Render rebuilds. To change the limits, set the command to
`python server.py --rate-runs 10 --global-runs 20` (Render: Settings > Docker Command).
