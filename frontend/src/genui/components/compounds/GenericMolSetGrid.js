import React from "react";
import { ResponsiveGrid } from '../../index';
import { Card } from 'reactstrap';

class GenericMolSetGrid extends React.Component {

  constructor(props) {
    super(props);

    this.cardComponent = this.props.cardComponent;
    this.newCardComponent = this.props.newCardComponent;
  }

  shouldComponentUpdate(nextProps, nextState, nextContext) {
    return this.props.molsets.length !== nextProps.molsets.length;
  }

  render() {
    const molsets = this.props.molsets;
    const headingText = this.props.headingText;

    const existing_cards = molsets.map(molset => ({
      id : molset.id,
      h : {"md" : 9, "sm" : 8},
      w : {"md" : 1, "sm" : 1},
      minH : {"md" : 3, "sm" : 3},
      data : molset
    }));
    const new_card = {
      id : "new-mol-set",
      h : {"md" : 7, "sm" : 6},
      w : {"md" : 1, "sm" : 1},
      minH : {"md" : 3, "sm" : 3},
      data : {}
    };

    const CardComponent = this.cardComponent;
    const NewCardComponent = this.newCardComponent;
    return (
      <React.Fragment>
        <h1>{headingText ? headingText : this.props.currentMolsetClass}</h1>
        <hr/>
        <ResponsiveGrid
          items={[new_card].concat(existing_cards)}
          rowHeight={75}
          mdCols={2}
          smCols={1}
          gridID={`${this.props.currentMolsetClass}-grid-layout`}
        >
          {
            [(
              <Card key={new_card.id} id={new_card.id}>
                <NewCardComponent {...this.props} handleCreateNew={this.props.handleAddMolSet}/>
              </Card>
            )].concat(existing_cards.map(
              item => (
                <Card key={item.id.toString()}>
                  <CardComponent
                    {...this.props}
                    molset={item.data}
                    onMolsetDelete={this.props.handleMolSetDelete}
                  />
                </Card>
              )
            ))
          }
        </ResponsiveGrid>
      </React.Fragment>
    )
  }
}

export default GenericMolSetGrid;