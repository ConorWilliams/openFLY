
#include "libfly/io/gsd.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"

void example_gsd() {
  //
  namespace fs = fly::system;

  fs::Box box(fly::Mat::Identity(), fly::Arr<bool>::Constant(true));

  fs::TypeMap<fly::Index> map(2);

  // Set up the map //

  fs::Supercell cell = fs::make_supercell<fly::Position>(box, map, 4);

  // Write some data to the supercell //

  fly::io::BinaryFile file("example.gsd", fly::io::create);

  file.commit([&] {
    file.write(cell.box());      //< Write the box to frame 0.
    file.write(cell.map());      //< Write the map to frame 0.
    file.write(cell.size());     //< Write the number of atoms to frame 0.
    file.write(fly::id_, cell);  //< Write the TypeID's of the atoms to frame 0.
    file.write(fly::r_, cell);   //< Write the positions of the atoms to frame 0.
  });

  // Do some data processing and compute new positions for each atom //

  file.commit([&] {
    file.write(fly::r_, cell);  //< Write the positions of the atoms to frame 1.

    // We do not need to duplicate box, map, etc. as they have not changed //
  });
}
