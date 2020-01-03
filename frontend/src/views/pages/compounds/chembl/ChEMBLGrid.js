import React, { Component } from 'react';
import { ResponsiveGrid } from '../../../../genui';
import { Card } from 'reactstrap';
import ChEMBLCard from './ChEMBLCard';
import ChEMBLCardNew from './ChEMBLCardNew';

class ChEMBLGrid extends Component {

  render() {
    const molsets = this.props.molsets;

    const existing_cards = molsets.map(molset => ({
      id : molset.id,
      h : {"md" : 6, "sm" : 6},
      w : {"md" : 1, "sm" : 1},
      minH : {"md" : 3, "sm" : 3},
      data : molset
    }));
    const new_card = {
      id : "new-mol-set",
      h : {"md" : 5, "sm" : 5},
      w : {"md" : 1, "sm" : 1},
      minH : {"md" : 3, "sm" : 3},
      data : {}
    };

    return (
      <React.Fragment>
        <h1>ChEMBL Compounds</h1>
        <hr/>
        <ResponsiveGrid
          items={existing_cards.concat(new_card)}
          rowHeight={75}
          mdCols={2}
          smCols={1}
        >
          {
            existing_cards.map(
              item => (
                <Card key={item.id.toString()}>
                  <ChEMBLCard {...this.props} molset={item.data}/>
                </Card>
              )
            ).concat([(
              <Card key={new_card.id} id={new_card.id}>
                <ChEMBLCardNew {...this.props}/>
              </Card>
            )])
          }
        </ResponsiveGrid>
      </React.Fragment>
    )
  }
}

export default ChEMBLGrid;