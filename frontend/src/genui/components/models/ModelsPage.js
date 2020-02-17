import React from 'react';
import { ComponentWithObjects, ModelGrid } from '../../index';

class ModelsPage extends React.Component {

  constructor(props) {
    super(props);

    this.headerComponent = this.props.headerComponent;

    this.state = {
      selectedToAdd : null,
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
    const selectedToAdd = this.state.selectedToAdd ? this.state.selectedToAdd : this.props.selectedToAdd;

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
                newModelComponent={this.state.newModelComponent}
                newCardSetup={this.state.newCardSetup}
                cardSetup={this.state.cardSetup}
                models={models[this.props.modelClass]}
                chosenAlgorithm={selectedToAdd}
                handleAddModel={
                  (...args) => {
                    this.setState({selectedToAdd : this.props.selectedToAdd});
                    return handleAddModel(...args)
                  }
                }
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