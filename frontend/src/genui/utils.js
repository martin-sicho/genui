
export function groupByMolset(mols, molsets, pushCallBack=(item, itemIdx) => item) {
  let molsetsIDs = molsets ? molsets.map(molset => molset.id) : null;
  const molsGroupedbyMolSet = {};
  mols.forEach((mol, index) => {
    // find the relevant molecule sets (only those on the map)
    const providers = mol.providers;
    const validProvidersIDs = [];
    providers.forEach((provider) => {
      if (molsetsIDs && molsetsIDs.includes(provider)) {
        validProvidersIDs.push(provider);
      } else {
        validProvidersIDs.push(provider);
      }
    });

    // assign the molecule to the correct molsets in the grouped representation
    validProvidersIDs.forEach(validProviderID => {
      if (!molsGroupedbyMolSet.hasOwnProperty(validProviderID)){
        molsGroupedbyMolSet[validProviderID] = []
      }
      molsGroupedbyMolSet[validProviderID].push(pushCallBack(mol, index));
    })
  });

  return molsGroupedbyMolSet;
}