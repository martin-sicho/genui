import React from "react";
import { GenericMolSetGrid } from '../../../../genui';
import GeneratedCard from './GeneratedCard';
import GeneratedCardNew from './GeneratedCardNew';

function GeneratedGrid(props) {
  const listUrl = new URL('generated/', props.apiUrls.compoundSetsRoot);
  return (
    <GenericMolSetGrid
      {...props}
      headingText="Generated Compounds"
      cardComponent={GeneratedCard}
      newCardComponent={GeneratedCardNew}
      molsetListUrl={listUrl}
    />
  )
}

export default GeneratedGrid