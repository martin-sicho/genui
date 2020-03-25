import { MolSetsTabs, MolsToMolSetGroups } from '../../../../genui';
import React from 'react';

export default function SelectedList(props) {
  const selectedMols = props.selectedMols;

  return (
    <MolsToMolSetGroups
      {...props}
      mols={selectedMols}
    >
      {
        (groups) => {
          return (
            <MolSetsTabs
              {...props}
              groupedMols={groups}
              mols={selectedMols}
            />
          );
        }
      }
    </MolsToMolSetGroups>
  );
}