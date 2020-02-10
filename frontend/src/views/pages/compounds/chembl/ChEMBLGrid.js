import React from 'react';
import {GenericMolSetGrid} from '../../../../genui';
import ChEMBLCard from './ChEMBLCard';
import ChEMBLCardNew from './ChEMBLCardNew';

class ChEMBLGrid extends React.Component {

  render() {
    return (
      <GenericMolSetGrid
        {...this.props}
        cardComponent={ChEMBLCard}
        newCardComponent={ChEMBLCardNew}
      />
    )
  }
}

export default ChEMBLGrid;