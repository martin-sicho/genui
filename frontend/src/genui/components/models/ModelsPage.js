import React from 'react';
import { ComponentWithObjects, ModelGrid } from '../../index';

class ModelsPage extends React.Component {

  constructor(props) {
    super(props);

    this.headerComponent = this.props.headerComponent;

    this.state = {
      selectedToAdd : this.props.selectedToAdd,
      newModelComponent : this.props.newModelComponent,
      newCardSetup : this.props.newCardSetup ? this.props.newCardSetup : { // TODO: it would be wiser to store all properties that will be passed to the card component rather than just this
        h : {"md" : 15, "sm" : 15},
        w : {"md" : 1, "sm" : 1},
        minH : {"md" : 3, "sm" : 3},
      },
      cardSetup : this.props.cardSetup ? this.props.cardSetup : { // TODO: it would be wiser to store all properties that will be passed to the card component rather than just this
        h : {"md" : 12, "sm" : 12},
        w : {"md" : 1, "sm" : 1},
        minH : {"md" : 3, "sm" : 3},
      }
    }
  }

  handleAddNew = (model, newModelComponent, cardSetup) => {
    this.setState((prevState) => {
      return {
        selectedToAdd : model,
        newModelComponent : newModelComponent ? newModelComponent : prevState.newModelComponent,
        newCardSetup : cardSetup ? cardSetup : prevState.newCardSetup,
      }
    })
  };

  componentDidMount() {
    const HeaderComp = this.headerComponent;
    if (HeaderComp) {
      this.props.onHeaderChange(<HeaderComp {...this.props} addChoices={this.props.algorithmChoices} onModelAdd={this.handleAddNew}/>);
    }
  }

  render() {
    const selectedToAdd = this.props.selectedToAdd ? this.props.selectedToAdd : this.state.selectedToAdd;
    const newModelComponent = this.props.newModelComponent ? this.props.newModelComponent : this.state.newModelComponent;
    const newCardSetup = this.props.newCardSetup ? this.props.newCardSetup : this.state.newCardSetup;
    const cardSetup = this.props.cardSetup ? this.props.cardSetup : this.state.cardSetup;

    return (
      <div className="models-page">
        <ComponentWithObjects
          {...this.props}
          emptyClassName={this.props.modelClass}
          objectListURL={this.props.listURL}
          render={
            (models, handleAddModelList, handleAddModel, handleModelDelete) => {
              return <ModelGrid
                {...this.props}
                chosenAlgorithm={selectedToAdd}
                newModelComponent={newModelComponent}
                newCardSetup={newCardSetup}
                cardSetup={cardSetup}
                models={models[this.props.modelClass]}
                handleAddModel={handleAddModel}
                handleModelDelete={handleModelDelete}
              />
            }
          }
        />
      </div>
    );
  }
}

export default ModelsPage;