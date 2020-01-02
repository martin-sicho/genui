import React, { Component } from 'react';
import { ResponsiveGrid } from '../../../../genui';
import { Card, CardHeader, CardBody, CardFooter, Button } from 'reactstrap';
import ChEMBLCard from './ChEMBLCard';

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
      h : {"md" : 3, "sm" : 3},
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
          rowHeight={100}
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
                <CardHeader>Card with Footer</CardHeader>
                <CardBody>
                  Lorem ipsum dolor sit amet, consectetur adipiscing elit. Morbi id neque quam. Aliquam-a-0 m-b-smitudin
                  egestas dui nec, fermentum diam. Vivamus vel tincidunt libero, vitae elementu
                </CardBody>
                <CardFooter>
                  <Button color="success">Add</Button> <Button>Cancel</Button>
                </CardFooter>
              </Card>
            )])
          }
        </ResponsiveGrid>
      </React.Fragment>
    )
  }
}

export default ChEMBLGrid;