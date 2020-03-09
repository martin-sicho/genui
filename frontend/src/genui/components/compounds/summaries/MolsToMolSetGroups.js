import { groupByMolset } from '../../../utils';

export default function MolsToMolSetGroups(props) {
  const mols = props.mols;
  const molsets = props.molsets;
  const groups = groupByMolset(mols, molsets, (mol, idx) => {
    return {
      mol: mol
    }
  });

  return props.children(groups);
}