

## Creatiing a new patch


### Topology
Create Topology using toposetDict

```bash
topoSet -dict topoSetDict.tank
```

Where the **topoSetDict.tank**  example creates a cylinder with axis defines between points p1 and p2 and radius 26

```cpp
actions
(
    {
        name    tank;
        type    faceSet;
        action  new;
        source  cylinderToFace;
        p1 (4736 2355 5.0);
        p2 (4736 2355 40.0);
        radius 26;
    }
);
```
This will select of cell faces in the domain that intersects with the cylinder and it will store it faceSet **tank**


### Restrict face

The faces selected before some of them are internal 

```bash
topoSet -dict topoSetDict.uniontank
```

```cpp
actions
(
  // create temp set with terrain faces (this is done as this mods temp1)
  {
    name    temp;
    type    faceSet;
    action  new;
    source  patchToFace;
    sourceInfo
      {
      name terrain;
      }
  }
  // keep in tank only faces intersecting with temp (terrain)    
  {
    name    tank;
    type    faceSet;
    action  subset;
    source  faceToFace;
          sourceInfo
          {
          set temp;
          option all;
          }
  }
);
```


### Create patch



```cpp
    {
        // Name of new patch
        name tankpatch;

        // Dictionary to construct new patch from
        patchInfo
        {
            type patch;
        }

        // How to construct: either from 'patches' or 'set'
        constructFrom set;

        // If constructFrom = set : name of faceSet
        set tank;
    }
```

```bash
createPatch -overwrite
```

