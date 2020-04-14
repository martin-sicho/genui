import React from "react";
import {Responsive, WidthProvider} from "react-grid-layout";

class ResponsiveGrid extends React.Component {

    getLayout = (items, cols, config) => {
        let row_id = 0;
        let col_id = 0;
        return items.map(item => {
            // console.log(item);
            const position = {
                i: item.id.toString()
                , x: col_id
                , y: row_id
                , w: item.w[config]
                , h: item.h[config]
                // , minW: 1
                // , maxW: 1
                , minH: item.minH[config]
            };
            if (col_id % cols) {
                row_id += 1;
                col_id = 0;
            } else {
                col_id += 1;
            }
            return position;
        });
    };


    render() {
        const items = this.props.items;
        const mdCols = (this.props.mdCols === undefined) ? 2 : this.props.mdCols;
        const smCols = (this.props.smCols === undefined) ? 1 : this.props.smCols;
        const mdBreak = (this.props.mdBreak === undefined) ? 992 : this.props.mdBreak;
        const smBreak = (this.props.smBreak === undefined) ? 480 : this.props.smBreak;

        const layouts = {
            md: this.getLayout(items, mdCols, "md")
            , sm: this.getLayout(items, smCols, "sm")
        };
        const ResponsiveGridLayout = WidthProvider(Responsive);
        // console.log(layouts);
        return (
            <ResponsiveGridLayout
                className={this.props.gridID}
                id={this.props.gridID}
                layouts={layouts}
                breakpoints={{md: mdBreak, sm: smBreak}}
                cols={{md: mdCols, sm: smCols}}
                rowHeight={this.props.rowHeight}
                autoSize={true}
                verticalCompact={true}
                draggableCancel='.unDraggable'
            >
                {this.props.children}
            </ResponsiveGridLayout>
        )
    }
}

export default ResponsiveGrid;