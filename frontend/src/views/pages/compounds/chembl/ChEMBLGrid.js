import React, { Component } from 'react';
import { ResponsiveGrid, TaskAwareComponent } from '../../../../genui';
import { Card } from 'reactstrap';
import ChEMBLCard from './ChEMBLCard';
import ChEMBLCardNew from './ChEMBLCardNew';

class ChEMBLGrid extends Component {

  constructor(props) {
    super(props);

    this.state = {
      molsets : this.props.molsets
    }
  }

  render() {
    const molsets = this.state.molsets;

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
                  <TaskAwareComponent
                    handleResponseErrors={this.props.handleResponseErrors}
                    tasksURL={new URL(`${item.data.id}/tasks/all/`, this.props.apiUrls.compoundSetsRoot)}
                    render={
                      (taskInfo, onTaskUpdate) => (
                        <ChEMBLCard
                          {...this.props}
                          {...taskInfo}
                          onTaskUpdate={onTaskUpdate}
                          molset={item.data}
                          onMolsetDelete={this.props.handleMolSetDelete}
                        />
                      )
                    }
                  />
                </Card>
              )
            ).concat([(
              <Card key={new_card.id} id={new_card.id}>
                <ChEMBLCardNew {...this.props} handleCreateNew={this.props.handleAddMolSet}/>
              </Card>
            )])
          }
        </ResponsiveGrid>
      </React.Fragment>
    )
  }
}

export default ChEMBLGrid;