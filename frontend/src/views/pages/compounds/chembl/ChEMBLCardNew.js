import React from 'react';
import { ChEMBLCreateForm } from './CreateForm';
import { CardHeader } from 'reactstrap';

class ChEMBLCardNew extends React.Component {

  createMolSet = (data) => {
    console.log(data)
  };

  render() {
    return (
      <React.Fragment>
        <CardHeader>Download New Compound Set from ChEMBL</CardHeader>
        <ChEMBLCreateForm handleCreate={this.createMolSet}/>
      </React.Fragment>
    )
  }

}

export default ChEMBLCardNew;