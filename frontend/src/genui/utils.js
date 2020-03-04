
export function filterProviders(mol, allowedProviders) {
  // find the relevant molecule sets (only those on the map)
  const providers = mol.providers;
  const validProviders = [];
  providers.forEach((provider) => {
    const molset = allowedProviders.find(x => x.id === provider);
    if (allowedProviders && molset) {
      validProviders.push(molset);
    }
  });

  return validProviders;
}

export function groupByMolset(mols, molsets, pushCallBack=(item, itemIdx) => item) {
  const molsGroupedbyMolSet = {};
  mols.forEach((mol, index) => {
    const validProviders = filterProviders(mol, molsets);

    // assign the molecule to the correct molsets in the grouped representation
    validProviders.forEach(validProvider => {
      const id = validProvider.id.toString();
      if (!molsGroupedbyMolSet.hasOwnProperty(id)){
        molsGroupedbyMolSet[id] = []
      }
      molsGroupedbyMolSet[id].push(pushCallBack(mol, index));
    })
  });
  return molsGroupedbyMolSet;
}