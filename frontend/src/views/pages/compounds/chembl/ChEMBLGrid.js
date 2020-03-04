import React from 'react';
import {GenericMolSetGrid} from '../../../../genui';
import ChEMBLCard from './ChEMBLCard';
import ChEMBLCardNew from './ChEMBLCardNew';

class ChEMBLGrid extends React.Component {

  render() {
    const listUrl = new URL('chembl/', this.props.apiUrls.compoundSetsRoot);
    return (
      <GenericMolSetGrid
        {...this.props}
        headingText="ChEMBL Compounds"
        cardComponent={ChEMBLCard}
        newCardComponent={ChEMBLCardNew}
        molsetListUrl={listUrl}
      />
    )
  }
}

export default ChEMBLGrid;