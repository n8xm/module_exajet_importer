# module_exajet_importer

To run in the example viewer
#Render the hexadron with particle

```bash
./ospExampleViewer --module exajet_import \
 --import:bin:/usr/sci/data/ospray/exajet-d12/hexas.bin
```



#Render the jet data with OSPRay unstructure mesh

```bash
./ospExampleViewer --module exajet_import \
  --import:jetunstr:<path to data>/hexas.bin \
  <path to data>/surfaces_vtp/*vtp
```

