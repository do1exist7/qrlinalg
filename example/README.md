# Complete lifecycle example

`example.f90` demonstrates a compact real-symmetric qrlinalg workflow:

1. allocate and initialize a state with spare capacity;
2. factorize and solve at one shift;
3. measure the action residual between physical matrices and stored factors;
4. replace, append, and delete symmetric rows and columns;
5. refactorize the current matrices when the shift changes;
6. deallocate the complete state;
7. allocate a fresh state with a different capacity and shift.

Build and run it from the repository root:

```sh
make example PREC=8
./build/release-wp8/qrlinalg_example
```

`PREC=10` and `PREC=16` select the other supported working kinds.

The state object is declared `allocatable` only to make explicit lifetime
control visible. Intrinsic Fortran `deallocate(state)` recursively releases
all allocatable private components owned by the state. A non-allocatable local
state receives the same automatic cleanup when its scope ends.

Calling `state%initialize(new_capacity, info)` is also safe without an explicit
deallocation. Initialization first releases the old factors and workspaces,
then builds an empty state for the new capacity. The shift is selected later
by `factorize_fresh`; initialization itself has no shift argument.

A shift change alone does not require reinitialization. Call
`factorize_fresh(H, S, new_shift, info)` on the existing initialized state.
Changing the shift is a dense change to `H-shift*S`, so replacement, append,
and deletion updates cannot be used for it.
