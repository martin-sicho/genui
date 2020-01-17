import React from "react";
import { ComponentWithMolSets } from '../../../genui';
// import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';

// function HeaderNav(props) {
//   return (<UncontrolledDropdown nav inNavbar>
//     <DropdownToggle nav caret>
//       Actions
//     </DropdownToggle>
//     <DropdownMenu right>
//       <UncontrolledDropdown>
//         <DropdownToggle nav>Add New...</DropdownToggle>
//         <DropdownMenu>
//           {
//             props.molSetChoices.map(choice =>
//               (<DropdownItem
//                 key={choice}
//                 onClick={() => {props.onMolSetChoice(choice, [])}}
//               >
//                 {choice}
//               </DropdownItem>)
//             )
//           }
//         </DropdownMenu>
//       </UncontrolledDropdown>
//     </DropdownMenu>
//   </UncontrolledDropdown>)
// }

class ModelsPage extends React.Component {

  componentDidMount() {
    // this.props.onHeaderChange(<HeaderNav {...this.props} molSetChoices={Object.keys(this.CLASS_TO_COMPONENT)} onMolSetChoice={this.props.handleAddMolSetList}/>);
  }

  render() {

    console.log(this.props.compoundSets);

    return (
      <div className="models-grid">
        This will be the grid with models
      </div>
    );
  }
}

function Models(props) {
  return (
    <ComponentWithMolSets
      {...props}
      render={
        (
          compoundSets
        ) => {
          return (<ModelsPage
            {...props}
            compoundSets={compoundSets}
          />)
        }
      }
    />
  )
}

export default Models;