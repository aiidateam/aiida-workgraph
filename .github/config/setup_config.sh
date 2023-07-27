verdi computer setup -n --config localhost-setup.yaml
verdi computer configure core.local localhost -n --config localhost-config.yaml
verdi code create core.code.installed -n --config code-pw.yaml --filepath-executable $(which pw.x)
verdi code create core.code.installed -n --config code-dos.yaml --filepath-executable $(which dos.x)
verdi code create core.code.installed -n --config code-projwfc.yaml --filepath-executable $(which projwfc.x)
verdi code create core.code.installed -n --config code-add.yaml --filepath-executable $(which bash)
aiida-pseudo install sssp -v 1.2 -x PBEsol -p efficiency
